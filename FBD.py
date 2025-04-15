import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle, Polygon, Arc
from ipywidgets import widgets
from IPython.display import display, clear_output

def draw_beam_structure(beam_type, length, load_type, load_value, load_position=None, load_span=None, moment_position=None):
    fig, ax = plt.subplots(figsize=(12, 4))

    # Draw the beam
    ax.plot([0, length], [0, 0], 'k-', lw=6)

    # Draw supports
    if beam_type == 'Cantilever':
        ax.add_patch(Rectangle((-0.5, -1), 0.5, 2, color='k'))
        ax.text(0, -1.5, 'Fixed', ha='center', fontsize=10)
    elif beam_type == 'Simply Supported':
        # Pin support (left)
        ax.plot(0, -0.5, 'ko', markersize=15)
        ax.plot([0, 0], [-0.5, -1.5], 'k-', lw=2)
        ax.add_patch(Polygon([[-0.5, -1.5], [0.5, -1.5], [0, -2.5]], color='k'))
        ax.text(0, -3, 'Pin', ha='center', fontsize=10)
        # Roller support (right)
        ax.plot(length, -0.5, 'ko', markersize=15)
        ax.plot([length, length], [-0.5, -1.5], 'k-', lw=2)
        ax.plot([length-0.5, length+0.5], [-2, -2], 'k-', lw=2)
        ax.add_patch(Polygon([[length-0.5, -1.5], [length+0.5, -1.5], [length, -2.5]], color='k'))
        ax.text(length, -3, 'Roller', ha='center', fontsize=10)
    elif beam_type == 'Fixed-Fixed':
        ax.add_patch(Rectangle((-0.5, -1), 0.5, 2, color='k'))
        ax.add_patch(Rectangle((length, -1), 0.5, 2, color='k'))
        ax.text(0, -1.5, 'Fixed', ha='center', fontsize=10)
        ax.text(length, -1.5, 'Fixed', ha='center', fontsize=10)

    # Draw loads
    if load_type == 'UDL':
        start, end = load_span if load_span else (0, length)
        num_arrows = int((end - start)/0.5) + 1
        for pos in np.linspace(start, end, num_arrows):
            ax.arrow(pos, 1.5, 0, -1.2, head_width=0.2, head_length=0.3, fc='r', ec='r', width=0.05)
        ax.plot([start, end], [1.5, 1.5], 'r-', lw=2)
        ax.text((start+end)/2, 1.8, f'{load_value} kN/m', ha='center', color='r')
    elif load_type == 'Point Load':
        ax.arrow(load_position, 2, 0, -1.6, head_width=0.3, head_length=0.4, fc='r', ec='r', width=0.1)
        ax.text(load_position, 2.2, f'{load_value} kN', ha='center', color='r')
    elif load_type == 'Moment':
        # Draw moment as a curved arrow
        arc = Arc((moment_position, 0), 1.5, 1.5, theta1=0, theta2=270, color='purple', lw=2)
        ax.add_patch(arc)
        # Arrowhead for the moment
        ax.arrow(moment_position+0.75, 0, -0.01, 0.01, head_width=0.15, head_length=0.15, fc='purple', ec='purple')
        ax.text(moment_position, 1.2, f'M={load_value} kNm', ha='center', color='purple', fontsize=12)

    ax.set_xlim(-1, length+1)
    ax.set_ylim(-4, 3)
    ax.axis('off')
    plt.title(f'{beam_type} Beam with {load_type}', fontsize=14)
    plt.tight_layout()
    plt.show()

# Widgets
beam_type_widget = widgets.Dropdown(
    options=['Cantilever', 'Simply Supported', 'Fixed-Fixed'],
    value='Cantilever',
    description='Beam Type:'
)
length_widget = widgets.FloatSlider(
    value=10.0, min=5.0, max=20.0, step=0.5, description='Length (m):'
)
load_type_widget = widgets.Dropdown(
    options=['UDL', 'Point Load', 'Moment'],
    value='UDL',
    description='Load Type:'
)
load_value_widget = widgets.FloatSlider(
    value=10.0, min=1.0, max=50.0, step=1.0, description='Load Value:'
)
load_position_widget = widgets.FloatSlider(
    value=5.0, min=0.0, max=20.0, step=0.5, description='Load Position (m):', layout={'visibility': 'hidden'}
)
udl_start_widget = widgets.FloatSlider(
    value=0.0, min=0.0, max=20.0, step=0.5, description='UDL Start (m):'
)
udl_end_widget = widgets.FloatSlider(
    value=10.0, min=0.0, max=20.0, step=0.5, description='UDL End (m):'
)
moment_value_widget = widgets.FloatSlider(
    value=10.0, min=-50.0, max=50.0, step=5.0, description='Moment Value (kNm):', layout={'visibility': 'hidden'}
)
moment_position_widget = widgets.FloatSlider(
    value=5.0, min=0.0, max=20.0, step=0.5, description='Moment Position (m):', layout={'visibility': 'hidden'}
)

def update_visibility(change):
    if change['new'] == 'Point Load':
        load_position_widget.layout.visibility = 'visible'
        udl_start_widget.layout.visibility = 'hidden'
        udl_end_widget.layout.visibility = 'hidden'
        moment_value_widget.layout.visibility = 'hidden'
        moment_position_widget.layout.visibility = 'hidden'
    elif change['new'] == 'UDL':
        load_position_widget.layout.visibility = 'hidden'
        udl_start_widget.layout.visibility = 'visible'
        udl_end_widget.layout.visibility = 'visible'
        moment_value_widget.layout.visibility = 'hidden'
        moment_position_widget.layout.visibility = 'hidden'
    else:  # Moment
        load_position_widget.layout.visibility = 'hidden'
        udl_start_widget.layout.visibility = 'hidden'
        udl_end_widget.layout.visibility = 'hidden'
        moment_value_widget.layout.visibility = 'visible'
        moment_position_widget.layout.visibility = 'visible'

def update_plot(*args):
    beam_type = beam_type_widget.value
    length = length_widget.value
    load_type = load_type_widget.value
    # Update slider limits
    for widget in [load_position_widget, udl_start_widget, udl_end_widget, moment_position_widget]:
        widget.max = length
        widget.value = min(widget.value, length)
    clear_output(wait=True)
    display(widgets.VBox([
        beam_type_widget,
        length_widget,
        load_type_widget,
        load_value_widget,
        moment_value_widget,
        load_position_widget,
        udl_start_widget,
        udl_end_widget,
        moment_position_widget
    ]))
    if load_type == 'Point Load':
        draw_beam_structure(beam_type, length, load_type, load_value_widget.value, load_position=load_position_widget.value)
    elif load_type == 'UDL':
        draw_beam_structure(beam_type, length, load_type, load_value_widget.value, load_span=(udl_start_widget.value, udl_end_widget.value))
    else:  # Moment
        draw_beam_structure(beam_type, length, load_type, moment_value_widget.value, moment_position=moment_position_widget.value)

load_type_widget.observe(update_visibility, names='value')
for widget in [beam_type_widget, length_widget, load_type_widget, load_value_widget, load_position_widget, udl_start_widget, udl_end_widget, moment_value_widget, moment_position_widget]:
    widget.observe(update_plot, names='value')

update_visibility({'new': load_type_widget.value})
display(widgets.VBox([
    beam_type_widget,
    length_widget,
    load_type_widget,
    load_value_widget,
    moment_value_widget,
    load_position_widget,
    udl_start_widget,
    udl_end_widget,
    moment_position_widget
]))
update_plot()
