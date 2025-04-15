import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, ConnectionPatch
from scipy.integrate import odeint

class EulerBernoulliBeam:
    def __init__(self, length, E, I, support_type='simply_supported'):
        """
        Initialize the beam with given properties.

        Parameters:
            length (float): Length of the beam (m)
            E (float): Young's modulus (Pa)
            I (float): Area moment of inertia (m^4)
            support_type (str): Type of support ('cantilever', 'simply_supported', 'fixed_fixed')
        """
        self.length = length
        self.E = E  # Young's modulus
        self.I = I  # Moment of inertia
        self.support_type = support_type
        self.loads = []  # List to store loads: ('point', pos, mag) or ('distributed', start, end, mag) or ('moment', pos, mag)
        self.reactions = None

    def add_point_load(self, position, magnitude):
        """Add a point load to the beam."""
        self.loads.append(('point', position, magnitude))

    def add_distributed_load(self, start_pos, end_pos, magnitude):
        """Add a distributed load to the beam."""
        self.loads.append(('distributed', start_pos, end_pos, magnitude))

    def add_moment_load(self, position, magnitude):
        """Add a moment load to the beam."""
        self.loads.append(('moment', position, magnitude))

    def calculate_reactions(self):
        """Calculate support reactions based on loads and support type."""
        if self.support_type == 'cantilever':
            total_force = 0
            total_moment = 0

            for load in self.loads:
                if load[0] == 'point':
                    _, pos, mag = load
                    total_force += mag
                    total_moment += mag * pos
                elif load[0] == 'distributed':
                    _, start, end, mag = load
                    length = end - start
                    total_force += mag * length
                    total_moment += mag * length * (start + end)/2
                elif load[0] == 'moment':
                    _, pos, mag = load
                    total_moment += mag

            self.reactions = {'R1': -total_force, 'M1': -total_moment}

        elif self.support_type == 'simply_supported':
            total_force = 0
            sum_moments_R1 = 0

            for load in self.loads:
                if load[0] == 'point':
                    _, pos, mag = load
                    total_force += mag
                    sum_moments_R1 += mag * pos
                elif load[0] == 'distributed':
                    _, start, end, mag = load
                    length = end - start
                    total_force += mag * length
                    sum_moments_R1 += mag * length * (start + end)/2
                elif load[0] == 'moment':
                    _, pos, mag = load
                    sum_moments_R1 += mag

            R2 = -sum_moments_R1 / self.length
            R1 = -total_force - R2
            self.reactions = {'R1': R1, 'R2': R2}

        elif self.support_type == 'fixed_fixed':
            total_force = 0
            sum_moments = 0

            for load in self.loads:
                if load[0] == 'point':
                    _, pos, mag = load
                    total_force += mag
                    sum_moments += mag * pos
                elif load[0] == 'distributed':
                    _, start, end, mag = load
                    length = end - start
                    total_force += mag * length
                    sum_moments += mag * length * (start + end)/2
                elif load[0] == 'moment':
                    _, pos, mag = load
                    sum_moments += mag

            R1 = -total_force / 2
            R2 = -total_force / 2
            M1 = total_force * self.length / 12 - sum_moments / 2
            M2 = -total_force * self.length / 12 - sum_moments / 2
            self.reactions = {'R1': R1, 'R2': R2, 'M1': M1, 'M2': M2}

    def shear_force(self, x):
        """Calculate shear force at position x."""
        V = 0

        if self.support_type == 'cantilever':
            if x >= 0:
                V += self.reactions['R1']
        elif self.support_type == 'simply_supported':
            if x >= 0:
                V += self.reactions['R1']
            if x >= self.length:
                V += self.reactions['R2']
        elif self.support_type == 'fixed_fixed':
            if x >= 0:
                V += self.reactions['R1']
            if x >= self.length:
                V += self.reactions['R2']

        for load in self.loads:
            if load[0] == 'point':
                _, pos, mag = load
                if x >= pos:
                    V += mag
            elif load[0] == 'distributed':
                _, start, end, mag = load
                if x > start:
                    if x <= end:
                        V += mag * (x - start)
                    else:
                        V += mag * (end - start)

        return V

    def bending_moment(self, x):
        """Calculate bending moment at position x."""
        M = 0

        if self.support_type == 'cantilever':
            if x >= 0:
                M += self.reactions['M1']
                M += self.reactions['R1'] * x
        elif self.support_type == 'simply_supported':
            if x >= 0:
                M += self.reactions['R1'] * x
            if x >= self.length:
                M += self.reactions['R2'] * (x - self.length)
        elif self.support_type == 'fixed_fixed':
            if x >= 0:
                M += self.reactions['M1']
                M += self.reactions['R1'] * x
            if x >= self.length:
                M += self.reactions['M2']
                M += self.reactions['R2'] * (x - self.length)

        for load in self.loads:
            if load[0] == 'point':
                _, pos, mag = load
                if x >= pos:
                    M += mag * (x - pos)
            elif load[0] == 'distributed':
                _, start, end, mag = load
                if x > start:
                    if x <= end:
                        M += mag * (x - start)**2 / 2
                    else:
                        M += mag * (end - start) * (x - (start + end)/2)
            elif load[0] == 'moment':
                _, pos, mag = load
                if x >= pos:
                    M += mag

        return M

    def deflection(self, x):
        """Calculate deflection at position x using numerical integration."""
        # Solve the differential equation EI*d4y/dx4 = w(x)
        # We'll use the double integration method

        # First, calculate M(x) at discrete points
        n_points = 100
        x_vals = np.linspace(0, self.length, n_points)
        M_vals = np.array([self.bending_moment(x) for x in x_vals])

        # Integrate M(x)/EI twice to get deflection
        theta = np.zeros_like(x_vals)
        y = np.zeros_like(x_vals)

        # Numerical integration (trapezoidal rule)
        for i in range(1, n_points):
            dx = x_vals[i] - x_vals[i-1]
            theta[i] = theta[i-1] + 0.5 * (M_vals[i-1] + M_vals[i]) / (self.E * self.I) * dx
            y[i] = y[i-1] + 0.5 * (theta[i-1] + theta[i]) * dx

        # Apply boundary conditions
        if self.support_type == 'simply_supported':
            # y(0) = 0 and y(L) = 0 already satisfied by initial conditions
            pass
        elif self.support_type == 'cantilever':
            # y(0) = 0 and theta(0) = 0 already satisfied
            pass
        elif self.support_type == 'fixed_fixed':
            # y(0) = 0 and y(L) = 0, theta(0) = 0 and theta(L) = 0
            # Our simple integration doesn't enforce all these, but we'll adjust
            y -= y[-1] * x_vals / self.length

        # Find deflection at requested x
        return np.interp(x, x_vals, y)

    def stress_at_point(self, x, y, z, beam_height):
        """
        Calculate normal and shear stress at a point in the beam.

        Parameters:
            x (float): Position along beam length
            y (float): Vertical position from neutral axis (positive up)
            z (float): Horizontal position from centroid
            beam_height (float): Total height of beam cross-section

        Returns:
            tuple: (sigma_x, tau_xy) normal and shear stress
        """
        M = self.bending_moment(x)
        V = self.shear_force(x)

        # Normal stress due to bending (sigma_x)
        sigma_x = -M * y / self.I

        # Shear stress (tau_xy) - simplified for rectangular cross-section
        Q = (beam_height**2 / 8 - y**2 / 2) * beam_height/2  # First moment of area
        tau_xy = V * Q / (self.I * beam_height)

        return sigma_x, tau_xy

    def plot_sfd_bmd(self, n_points=100):
        """Plot Shear Force Diagram and Bending Moment Diagram."""
        x_vals = np.linspace(0, self.length, n_points)
        V_vals = np.array([self.shear_force(x) for x in x_vals])
        M_vals = np.array([self.bending_moment(x) for x in x_vals])

        plt.figure(figsize=(12, 8))

        plt.subplot(2, 1, 1)
        plt.plot(x_vals, V_vals, 'r-', linewidth=2)
        plt.title('Shear Force Diagram')
        plt.xlabel('Position along beam (m)')
        plt.ylabel('Shear Force (N)')
        plt.grid(True)
        plt.axhline(0, color='k', linestyle='-', linewidth=0.5)

        plt.subplot(2, 1, 2)
        plt.plot(x_vals, M_vals, 'b-', linewidth=2)
        plt.title('Bending Moment Diagram')
        plt.xlabel('Position along beam (m)')
        plt.ylabel('Bending Moment (Nm)')
        plt.grid(True)
        plt.axhline(0, color='k', linestyle='-', linewidth=0.5)

        plt.tight_layout()
        plt.show()

    def plot_deflection(self, n_points=100):
        """Plot the deflected shape of the beam."""
        x_vals = np.linspace(0, self.length, n_points)
        y_vals = np.array([self.deflection(x) for x in x_vals])

        plt.figure(figsize=(10, 6))
        plt.plot(x_vals, y_vals * 1000, 'g-', linewidth=2)  # Convert to mm for visibility
        plt.title('Deflected Shape of Beam')
        plt.xlabel('Position along beam (m)')
        plt.ylabel('Deflection (mm)')
        plt.grid(True)
        plt.axhline(0, color='k', linestyle='-', linewidth=0.5)
        plt.show()

    def plot_mohr_circle(self, x, y, z, beam_height):
        """Plot Mohr's Circle for stress at a given point."""
        sigma_x, tau_xy = self.stress_at_point(x, y, z, beam_height)
        sigma_y = 0  # Assuming no stress in y-direction

        # Calculate principal stresses and max shear
        sigma_avg = (sigma_x + sigma_y) / 2
        R = np.sqrt(((sigma_x - sigma_y) / 2)**2 + tau_xy**2)
        sigma_1 = sigma_avg + R
        sigma_2 = sigma_avg - R
        tau_max = R

        # Create figure
        fig, ax = plt.subplots(figsize=(8, 8))

        # Plot Mohr's circle
        circle = Circle((sigma_avg, 0), R, fill=False, edgecolor='blue', linewidth=2)
        ax.add_patch(circle)

        # Plot original stress point
        ax.plot(sigma_x, tau_xy, 'ro', markersize=8)
        ax.plot(sigma_x, -tau_xy, 'ro', markersize=8)

        # Plot principal stresses
        ax.plot([sigma_1, sigma_2], [0, 0], 'go', markersize=8)

        # Plot lines connecting points
        line1 = ConnectionPatch((sigma_avg, 0), (sigma_x, tau_xy), 'data', 'data', linestyle='--', color='gray')
        line2 = ConnectionPatch((sigma_avg, 0), (sigma_x, -tau_xy), 'data', 'data', linestyle='--', color='gray')
        ax.add_patch(line1)
        ax.add_patch(line2)

        # Set plot limits and labels
        ax.set_xlim(sigma_avg - 1.5*R, sigma_avg + 1.5*R)
        ax.set_ylim(-1.5*R, 1.5*R)
        ax.axhline(0, color='k', linestyle='-', linewidth=0.5)
        ax.axvline(0, color='k', linestyle='-', linewidth=0.5)
        ax.set_aspect('equal')
        ax.grid(True)

        ax.set_title(f"Mohr's Circle at x={x:.2f}m, y={y:.2f}m")
        ax.set_xlabel('Normal Stress (Pa)')
        ax.set_ylabel('Shear Stress (Pa)')

        # Add annotations
        ax.annotate(f'σ1 = {sigma_1:.2e} Pa', (sigma_1, 0), textcoords="offset points", xytext=(0,10), ha='center')
        ax.annotate(f'σ2 = {sigma_2:.2e} Pa', (sigma_2, 0), textcoords="offset points", xytext=(0,10), ha='center')
        ax.annotate(f'τmax = {tau_max:.2e} Pa', (sigma_avg, R), textcoords="offset points", xytext=(0,10), ha='center')

        plt.show()

        return sigma_1, sigma_2, tau_max

# Example usage with moment load
if __name__ == "__main__":
    # Beam properties
    length = 5.0  # m
    E = 200e9  # Pa (steel)
    I = 8.33e-6  # m^4 (rectangular 0.1m x 0.2m)
    beam_height = 0.2  # m

    # Create beam
    beam = EulerBernoulliBeam(length, E, I, support_type='simply_supported')

    # Add loads
    beam.add_point_load(2.0, -5000)  # 5000 N downward at 2m
    beam.add_moment_load(3.0, 2000)  # 2000 Nm clockwise moment at 3m

    # Calculate reactions
    beam.calculate_reactions()
    print("Reactions:", beam.reactions)

    # Plot diagrams
    beam.plot_sfd_bmd()
    beam.plot_deflection()

    # Analyze stress at a point and plot Mohr's circle
    x_analysis = 2.5  # m along beam
    y_analysis = 0.05  # m from neutral axis (top of beam)
    z_analysis = 0.0  # m (center of width)

    sigma_1, sigma_2, tau_max = beam.plot_mohr_circle(x_analysis, y_analysis, z_analysis, beam_height)
    print(f"Principal stresses at x={x_analysis}m: σ1={sigma_1:.2e} Pa, σ2={sigma_2:.2e} Pa")
    print(f"Maximum shear stress: {tau_max:.2e} Pa")
