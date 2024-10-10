import dash
from dash import dcc, html
from dash.dependencies import Input, Output, MATCH, ALL
import dash.exceptions
import dash_bootstrap_components as dbc
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
matplotlib.use('Agg')
import io
import base64
import numpy as np


#app = dash.Dash(__name__, suppress_callback_exceptions=True)
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP],suppress_callback_exceptions=True)
server = app.server  # For deployment if needed


# Helper function to create the Matplotlib figure
def create_beam_figure(n_supports,spacings,left,right,location):

    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(14, 6))


    #spacing=5 #this is initial spacing used between all supports. will update later will callback

    # Beam length calculation (supports are spaced 5 ft apart initailly)
    beam_length = sum(spacings)

    # Draw the beam as a horizontal line
    ax[0].plot([0, beam_length], [0, 0], color='black', linewidth=5, label='Beam')

    #The supports are all going to be simple supports except potentially the end supports.
    support_types= ['Simple','Moment','None (Overhang)']
    #create list of support distances from beginning of beam:







    # Add support patches (triangle markers)
    current_distance=0
    for i,spacing in enumerate(spacings):
        current_distance += spacing
        x_position = current_distance
        if x_position==0: #far left support
            if left==support_types[1]:
               support=patches.Rectangle((-.06,-.25),.05,.5, edgecolor='black', facecolor='red')
            elif left==support_types[2]:
                support= patches.Rectangle((-.06,-.25),.05,.5, edgecolor='black', facecolor='red',visible=False)
            else:
                support = patches.Polygon(np.array(
                    [[x_position, -0.05], [x_position - .25, -0.25],
                     [x_position + .25, -0.25]]),
                    edgecolor='black', facecolor='red')
        elif x_position==beam_length: #far left support
            if right==support_types[1]:
               support=patches.Rectangle((x_position+.01,-.25),.05,.5, edgecolor='black', facecolor='red')
            elif right==support_types[2]:
                support= patches.Rectangle((x_position+.01,-.25),.05,.5, edgecolor='black', facecolor='red',visible=False)
            else:
                support = patches.Polygon(np.array(
                    [[x_position, -0.05], [x_position - .25, -0.25],
                     [x_position + .25, -0.25]]),
                    edgecolor='black', facecolor='red')
        else:
            support = patches.Polygon(np.array(
                [[x_position, -0.05], [x_position - .25, -0.25],
                 [x_position + .25, -0.25]]),
                edgecolor='black', facecolor='red')
        ax[0].add_patch(support)
        ax[0].text(x_position, -0.4, f'{chr(65 + i)}', ha='center', va='top')
        if not x_position==0:

            ax[0].annotate(
                '', xy=(x_position, -.75), xytext=(x_position-spacing, -.75),  # Starting and ending points of the arrow
                arrowprops=dict(arrowstyle='<->', lw=2)  # Two-way arrow style with linewidth
            )

            # Add text in the middle of the arrow
            ax[0].text(x_position-spacing/2, -.6, f'{spacing} ft')

    # Set the axis limits
    ax[0].set_xlim(-1, beam_length + 1)
    ax[0].set_ylim(-1, 1)


    # Set x-axis labels and title
    ax[0].set_xlabel('Position along beam (ft)')
    ax[0].set_title(f'Beam with {n_supports} Supports')


    #create influence line subplot

    x = np.arange(0, beam_length, .01)

    ax[1].plot(x, y)
    ax[1].set_xlim(-1, beam_length + 1)
    ax[1].set_ylim(-1, max(y))
    # Set x-axis labels and title
    ax[1].set_xlabel('unit load position on beam')
    ax[1].set_title(f'Influence @ x= {n_supports} ft') #will change this once slider value is incorporated

    # Adjust layout to prevent overlap
    plt.tight_layout()

    # Save the figure to a buffer
    buf = io.BytesIO()
    plt.savefig(buf, format='png', bbox_inches='tight')
    plt.close(fig)
    buf.seek(0)

    # Encode the image to base64 for rendering in Dash
    image_base64 = base64.b64encode(buf.read()).decode('utf-8')
    return "data:image/png;base64,{}".format(image_base64)

    #create influence line subplot
'''
def create_influence_chart(n_supports):
    fig, ax1 = plt.subplots(figsize=(15, 2))
    spacing = 5
    beam_length = (n_supports - 1) * spacing
    x=np.arange(0,beam_length,.01)
    y=np.array([i**2 for i in x])
    ax1.plot(x,y)

    # Save the figure to a buffer
    buf = io.BytesIO()
    plt.savefig(buf, format='png', bbox_inches='tight')
    plt.close(fig)
    buf.seek(0)

    # Encode the image to base64
    image_base64 = base64.b64encode(buf.read()).decode('utf-8')
    return "data:image/png;base64,{}".format(image_base64)

'''

# Layout for the Dash app
app.layout = html.Div([
    html.Div([html.H3('Beam Analysis')],style={'display': 'flex', 'align-items': 'center','justify-content': 'center', 'margin-bottom': '20px'}),
    html.Div([
        html.Label('Select number of supports:', style={'margin-right': '10px'}),
        dcc.Dropdown(
            id='support-dropdown',
            options=[{'label': str(i), 'value': i} for i in range(2, 11)],
            value=2,  # Default value
            clearable=False,
            style={'width': '120px'}
        )
    ], style={'display': 'flex', 'align-items': 'center','justify-content': 'center', 'margin-bottom': '20px'}),


# Container for dynamically generated input fields for custom spacing
    html.Div(id='input-boxes', style={'margin-bottom': '20px'}),

html.Br(),
html.Br(),

# Dropdown for far left support type
html.Div(
    children=[
        html.Div(
            children=[
                # Left-aligned label and dropdown
                html.Div(
                    children=[
                        html.Label("Left End Support Type:", style={'margin-right': '10px'}),
                        dcc.Dropdown(
                            id='left-dropdown',
                            options=['Simple', 'Moment', 'None (Overhang)'

                                     ],
                            value='Simple',
                            style={'width': '200px'}
                        ),
                    ],
                    style={'display': 'flex', 'align-items': 'center'}  # Flexbox for left items
                ),
                # Right-aligned label and dropdown
                html.Div(
                    children=[
                        html.Label("Right End Support Type:", style={'margin-right': '10px'}),
                        dcc.Dropdown(
                            id='right-dropdown',
                            options=['Simple', 'Moment', 'None (Overhang)'

                                     ],
                            value='Simple',
                            style={'width': '200px'}
                        ),
                    ],
                    style={'display': 'flex', 'align-items': 'center'}  # Flexbox for right items
                ),
            ],
            style={
                'display': 'flex',  # Use flexbox for layout
                'justify-content': 'space-between',  # Spread items to the left and right
                'align-items': 'center',  # Center items vertically
                'width': '100%',  # Full width container
                'padding': '10px'
            }
        ),
    ]
),

    # Image for the Matplotlib plot
    html.Img(id='beam-plot', style={'margin-top': '20px','width': '100%'}),

# Display the sum of spacings. this is really just a placeholder to confirm the spacing values get updated
    html.Div(id='output-sum', style={'font-size': '24px', 'font-weight': 'bold', 'margin-top': '20px'}),
# Display the influence position. this is really just a placeholder to confirm the slider value gets updated
    html.Div(id='influence-position', style={'font-size': '24px', 'font-weight': 'bold', 'margin-top': '10px'}),
# Slider to select influence position
    html.Br(),
    html.Br(),

            html.Label("Select the influence position:"),
            dcc.Slider(
                id='slider-x',
                min=0,
                max=10,  # Will be updated based on L
                step=.01,
                value=0,
                marks={i: str(i) for i in range(11)},
                updatemode='drag'  # Update value while dragging

            )

])

@app.callback(
    Output('input-boxes','children'),
    Input('support-dropdown', 'value'),

        )

def update_spacing_inputs(n_supports):
    inputs = []
    # Create input fields for the spacing between each pair of supports
    for i in range(n_supports - 1):
        inputs.append(html.Div([
            html.Label(f'Spacing between support {chr(65 + i)} and {chr(65 + i + 1)}:', style={'margin-right': '10px'}),
            dcc.Input(
                id={'type': 'input-box', 'index': i + 1},
                type='number',
                value=5,  # Default value is 5 ft
                min=0,
                style={'width': '100px', 'margin-right': '20px'}
            )
        ], style={'display': 'flex', 'align-items': 'center', 'justify-content': 'center'}))
    return inputs
# Callback to update the Matplotlib graph when the dropdown value changes
'''
@app.callback(
    Output('beam-plot', 'src',allow_duplicate=True),
    Input('support-dropdown', 'value'),
    prevent_initial_call=True
)
def update_beam_plot(n_supports):
    # Create the beam plot based on the number of supports
    return create_beam_figure(n_supports)
'''
@app.callback(
    Output('influence-position', 'children'),
    Input('slider-x', 'value')
)
def update_sum(pos):
    # get the value from the slider and display it
    Influence = pos
    return f'Influence Position: {Influence}'


# Callback to update the sum of the input values based on the number of input boxes generated. this is not uselful for program, just starting template for next callback

@app.callback(
    Output('output-sum', 'children'),
    Input({'type': 'input-box', 'index': ALL}, 'value')  # Dynamically match all input boxes
)
def update_sum(values):
    # Calculate the sum of the input values, treating None values as 0
    total_sum = sum(value if value is not None else 0 for value in values)
    return f'Total Sum: {total_sum}'

@app.callback(
    Output('beam-plot', 'src'),
    [Input('support-dropdown', 'value'),
     Input({'type': 'input-box', 'index': ALL}, 'value'),
     Input('left-dropdown','value'),
     Input('right-dropdown','value'),
     Input('slider-x', 'value')],
)
def update_beam(n_supports,values,left,right,location):
    # create a list of current support spacings

    spacings = [i if i is not None else 0 for i in values]
    spacings.insert(0,0)
    return create_beam_figure(n_supports,spacings,left,right,location)






# Run the Dash app
if __name__ == '__main__':
    app.run_server(debug=True)

