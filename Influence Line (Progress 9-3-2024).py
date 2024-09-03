import tkinter as tk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
import numpy as np
import mplcursors


class BeamInfluenceLine:
    def __init__(self, root):
        self.root = root
        self.root.title("Influence Line of a Beam")

        # Reset Button
        self.reset_button = tk.Button(root, text="Reset", command=self.reset)
        self.reset_button.place(x=.9,y=.5)

        # Create input field for number of supports
        self.supports_label = tk.Label(root, text="step 1: click on dropdown menu below to select number of supports. then press 'next'.")
        self.supports_label.pack()
        #
        def create_dropdown(options, master):
            variable = tk.StringVar(master)
            variable.set(options[0])  # default value
            dropdown = tk.OptionMenu(master, variable, *options)
            return dropdown, variable

        self.options = [2, 3, 4, 5]
        self.dropdown, self.variable = create_dropdown(self.options, root)

        self.dropdown.pack()


        # Create button to create spacing input fields
        self.create_spacing_button = tk.Button(root, text="Next", command=self.create_spacing_fields)
        self.create_spacing_button.pack()

        # Create plot
        self.figure = plt.Figure(figsize=(10, 5), dpi=100)
        self.axes1 = self.figure.add_subplot(211)
        self.axes2 = self.figure.add_subplot(212)
        self.axes1.set_title("Beam Diagram")
        self.axes2.set_title("Influence Line")
        self.axes2.set_xlabel("Beam Length")
        self.axes2.set_ylabel("Influence Line Value")
        self.figure.tight_layout(pad=2)

        #create a support option for the start (1) and end (2) of beam. also create statuses for on_click fucntion to access:
        self.bubble_S1 = (plt.Rectangle((.05, .15), 0.01,height=.05, facecolor='blue', edgecolor='black', transform=self.axes1.transAxes, visible=False))
        self.axes1.add_patch(self.bubble_S1)
        self.bubble_S1_status = True
        self.bubble_M1 = (plt.Rectangle((.05, .05), 0.01,height=.05, facecolor='white', edgecolor='black', transform=self.axes1.transAxes, visible=False))
        self.axes1.add_patch(self.bubble_M1)
        self.bubble_M1_status = False
        self.bubble_S2 = (plt.Rectangle((.975, .15), 0.01, height=.05, facecolor='blue', edgecolor='black',
                                       transform=self.axes1.transAxes, visible=False))
        self.bubble_S2_status = True
        self.axes1.add_patch(self.bubble_S2)
        self.bubble_M2 = (plt.Rectangle((.975, .05), 0.01, height=.05, facecolor='white', edgecolor='black',
                                        transform=self.axes1.transAxes, visible=False))
        self.axes1.add_patch(self.bubble_M2)
        self.bubble_M2_status = False

        self.canvas = FigureCanvasTkAgg(self.figure, master=root)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack()


        self.cursor = None

        # Create a label to display the coordinates
        coord_label = tk.Label(root, text="")
        coord_label.pack()

    # Function to update the coordinates label
        def update_coords(event):
            coord_label.config(text=f"({event.xdata}, {event.ydata})")

        # Bind the motion_notify_event to the update_coords function
        self.canvas.mpl_connect('motion_notify_event', update_coords)

    def create_spacing_fields(self):
        # Get number of supports
        self.info = tk.Label(root, text="step 2: determine the support types in the beam diagram, and the spacing between supports below. then click 'update plot'. ")
        self.info.pack()

        num_supports = self.variable.get()

        num_supports = int(num_supports)
        # Create input fields for spacing
        self.spacing_entries = []
        # create support names
        self.supports_label = ['A', 'B', 'C', 'D', 'E']

        # Create initial beam diagram for visual purposes
        beam = patches.Rectangle((0, 0.5), (num_supports - 1) * 5, 2, edgecolor='black', facecolor='grey')
        self.axes1.add_patch(beam)
        #add in supports
        for i in range(0, (num_supports)):
            if i >= 0:  # Check if spacings list is not empty
                support = patches.Polygon(np.array(
                    [[i * 5, 0.5], [i * 5 - .25, -0.5],
                     [i * 5 + .25, -0.5]]),
                    edgecolor='black', facecolor='red')
                self.axes1.add_patch(support)
                self.axes1.text(i * 5, -1, s=f'{self.supports_label[i]}')
        self.axes1.set_xlim(-1, (num_supports - 1) * 5 + 1)
        self.axes1.set_ylim(-3, 4)


        # make selection bubbles visible AND ADD IN TEXT
        self.axes1.text(.025, .15, 'S',color='blue', transform=self.axes1.transAxes)
        self.bubble_S1.set_visible(True)
        self.axes1.text(.025, .05, 'M',color='blue', transform=self.axes1.transAxes)
        self.bubble_M1.set_visible(True)
        self.axes1.text(.95, .15, 'S', color='blue', transform=self.axes1.transAxes)
        self.bubble_S2.set_visible(True)
        self.axes1.text(.95, .05, 'M', color='blue', transform=self.axes1.transAxes)
        self.bubble_M2.set_visible(True)

        self.canvas.draw()
        self.canvas.mpl_connect('button_press_event', self.on_click)

        for i in range(num_supports - 1):
            spacing_label = tk.Label(self.root,
                                     text=f"Spacing {self.supports_label[i]} - {self.supports_label[i + 1]} (ft):")
            spacing_label.pack()
            spacing_entry = tk.Entry(self.root, width=5, )
            spacing_entry.insert(0, 5)


            spacing_entry.pack()
            self.spacing_entries.append(spacing_entry)

        # Create button to update plot
        self.update_button = tk.Button(self.root, text="Update Plot", command=self.update_plot)
        self.update_button.pack()
        # Create Slider for influence line section of interest
        self.space = tk.Label(self.root, text="")
        self.space.pack()
        self.slider_label = tk.Label(self.root, text="step 3: drag slider below to the influence position on the beam, indicated by the dashed red line in the Beam Diagram. you can hover over the influence Line to see the influence value.")
        self.slider_label.pack()

        self.location_slider = tk.Scale(self.root, from_=0,
                                        to=1,
                                        resolution=0.01,
                                        orient=tk.HORIZONTAL,
                                        command=self.update_plot, length=900, sliderlength=10, tickinterval=.1)
        self.location_slider.pack()




    # Function to handle support click
    def on_click(self,event):
        if event.xdata is not None and event.ydata is not None:
            if -.4 < event.xdata < -.26 and -2 < event.ydata < -1.5:
                # Select option S
                self.bubble_S1.set_facecolor('blue')
                self.bubble_M1.set_facecolor('white')
                self.bubble_S1_status = True
                self.bubble_M1_status = False
            elif -.4 < event.xdata < -.26 and -2.7 < event.ydata < -2.2:
                # Select option M
                self.bubble_S1.set_facecolor('white')
                self.bubble_M1.set_facecolor('blue')
                self.bubble_S1_status = False
                self.bubble_M1_status = True
            elif 10.68 < event.xdata < 10.82 and -2 < event.ydata < -1.5:
                # Select option S
                self.bubble_S2.set_facecolor('blue')
                self.bubble_M2.set_facecolor('white')
                self.bubble_S2_status = True
                self.bubble_M2_status = False
            elif 10.68 < event.xdata < 10.82 and -2.7 < event.ydata < -2.2:
                # Select option M
                self.bubble_S2.set_facecolor('white')
                self.bubble_M2.set_facecolor('blue')
                self.bubble_S2_status = False
                self.bubble_M2_status = True


            self.canvas.draw()


    def update_plot(self, val=None):
            # Get input values
            num_supports = int(self.variable.get())
            spacings = [float(entry.get()) for entry in self.spacing_entries]
            support_positions = [sum(spacings[:i]) for i in range(num_supports + 1)]
            location = self.location_slider.get()

            # Clear previous plot
            self.axes1.clear()
            self.axes2.clear()
            self.axes1.set_title("Beam Diagram")
            self.axes2.set_title(f'Influence Line @x={round(location * sum(spacings), 2)}')
            self.axes1.set_xlabel("Beam Length")
            self.axes2.set_xlabel("Beam Length")
            self.axes2.set_ylabel("Shear influence")

            # Plot beam diagram
            beam = patches.Rectangle((0, 0.5), sum(spacings), 2, edgecolor='black', facecolor='grey')
            self.axes1.add_patch(beam)
            #arrow=self.axes1.arrow(0, 4, 0, -1, width=.1)
            for i, pos in enumerate(support_positions):
                if i <= len(spacings):  # Check if spacings list is not empty
                    if i == 0 and self.bubble_S1_status == False:
                        support = patches.Rectangle((-.150, -.5), .15, 4, edgecolor='black', facecolor='red')
                    elif i == len(spacings) and self.bubble_S2_status == False:
                        support = patches.Rectangle((pos, -.5), .15, 4, edgecolor='black', facecolor='red')
                    else:
                        support = patches.Polygon(np.array([[pos, 0.5], [pos - .25, -0.5], [pos + .25, -0.5]]),
                                                  edgecolor='black', facecolor='red')
                    self.axes1.add_patch(support)
                    self.axes1.text(pos, -1, s=f'{self.supports_label[i]}')

            self.axes1.axvline(x=location * sum(spacings), color='r', linestyle='--')  # Add vertical dashed red line
            self.axes1.set_xlim(-1, sum(spacings) + 1)
            self.axes1.set_ylim(-3, 4)
            # self.axes1.set_aspect('equal', adjustable='box')

            # Plot influence line

            if num_supports == 2:
                x = np.arange(0, sum(spacings), .01)
                y = [-i / sum(spacings) if i < location * sum(spacings) else (sum(spacings) - i) / sum(spacings) for i in x]
                self.axes2.text(x=location * sum(spacings), y=-(location * 1.2), s=f'{str(round(-location, 2))} (max.)')
                self.axes2.text(x=location * sum(spacings), y=(1 - location) * 1.2,
                                s=f'{str(round(1 - location, 2))} (max.)')

            if num_supports == 3:
                x = np.arange(0, sum(spacings), .01)
                # a beam with 3 or more supports is considered statically determinant, thus making the classic equilibrium
                # equations insufficient for any shear calculations. rather, before any influence line equations can be derived,
                # we must first solve a reaction in terms of the load location "x" and then can continue
                # with our new adjusted statically "determinant" structure

                # solve the middle support in terms of load location
                L = sum(spacings)  # total length of beam
                c = L * location  # location of interest for influence line diagram, determined by the slider
                L1 = spacings[0]  # distance of middle support from left end
                L2 = spacings[1]  # distance of middle support from right end

                # we will solve the reaction(x) of the middle support by superposition (removing the support
                # and solving the deflection at the middle support location, then solving the required P load at
                # the middle support to produce the same upwards deflection. the following list determines the
                # deflection of the middle support when the load is applied at every position "a" on the beam

                delta_Bx = np.array(
                    [(L - i) * L1 / (6 * L) * (L ** 2 - (L - i) ** 2 - L1 ** 2) if i >= L1 else i * L2 / (6 * L) * (
                            L ** 2 - i ** 2 - L2 ** 2) for i in x])
                # this is deflection of middle support as a function of the load location

                # delta_BB= P(L1^2)(L2^2)/3L , this is the deflection @ middle support when upward P load applied at the middle support
                # equating both equations allows us to solve for P (shear @influence location when load placed there):
                B_y = np.array([i / ((L1 ** 2) * (L2 ** 2) / (3 * L)) for i in delta_Bx])

                # Now that the middle reaction is determined, the other supports can be solved as functions of load position "a" as well.

                A_y = np.array([(-i * L2 + (L - j)) / L for i, j in zip(B_y, x)])
                C_y = np.array([(-i * L1 + j) / L for i, j in zip(B_y, x)])

                # so far we have created lists for the reactions of A, B , and C when the load is placed on the beam
                # we will use these reaction values to determine the positive and negative influence value at the slider location:
                # using simple statics:

                # I_loc_neg=np.array([location,A_y[np.where(x == location)[0][0]]-1])
                # I_loc_pos=np.array([location,A_y[np.where(x == location)[0][0]]])

                # now we can add the remaining known Influence values, which is when the load is placed at the supports:
                # if location == 0:
                #     I_Ay=np.array([0,1])
                #     I_By=np.array([spacings[0],0])
                #     I_Cy = np.array([spacings[0]+spacings[1],0])
                # if location==spacings[0]:
                #     I_Ay = np.array([0, 0])
                #     I_By = np.array([spacings[0],1])
                #     I_Cy = np.array([spacings[0] + spacings[1], 0])
                #  if location==spacings[0]+spacings[1]:
                #     I_Ay = np.array([0, 0])
                #     I_By = np.array([spacings[0],0])
                #     I_Cy = np.array([spacings[0] + spacings[1], 1])
                #
                # else:
                #     I_Ay = np.array([0, 0])
                #     I_By = np.array([spacings[0], 0])
                #     I_Cy = np.array([spacings[0] + spacings[1], 0])

                # Now that we have the Ay and By reactions for all load positions, we can determine the influence value
                # for "location"
                if c == 0:
                    y = np.array([i for i in A_y])
                elif c < L1:
                    y = np.array([i - 1 if j < c else i for i, j in zip(A_y, x)])
                elif c == L1:
                    y = np.array([i for i in B_y])
                elif c<L:
                    y = np.array([(i + k - 1) if j < c else (i + k) for i, k, j in zip(A_y, B_y, x)])
                else:
                    y = np.array([-i for i in C_y])

            # PLACEHOLDER FOR SUPPORTS 4 AND 5
            if num_supports > 3:
                x = np.arange(0, sum(spacings), .01)
                y = [1 for i in x]

            # Capture the line object
            line = self.axes2.plot(x, y)

            self.axes2.plot(x, y)
            self.axes2.axhline(y=0, color='r', linestyle='-')  # Add horizontal red line for x-axis
            self.axes2.set_xlim(-1, sum(spacings) + 1)
            self.axes2.set_ylim(-1.5, 1.5)

            #determine the load location for the point load display

            arrow= patches.FancyArrowPatch((0.5, 2.5), (.5, 3.5), arrowstyle='<-',mutation_scale=10, color='red', visible=False)
            self.axes1.add_patch(arrow)
            arrow_text=self.axes1.text(.5, 3.4, s='P=1')


            last_hover_x = 0.5

            def move_arrow(sel):
                global last_hover_x
                x, y = sel.target
                last_hover_x = x
                arrow.set_positions((last_hover_x, 2.5), (last_hover_x, 3.5))
                arrow.set_visible(True)
                arrow_text.set_visible(True)
                arrow_text.set_position((last_hover_x, 3.4))
                sel.annotation.set_text(f"x: {x:.2f}, y: {y:.2f}")
                self.figure.canvas.draw_idle()


            def remove_arrow(sel):
                arrow.set_xy = (last_hover_x,4)
                self.figure.canvas.draw_idle()





            # Clear old hovers
            if self.cursor:
                self.cursor.remove()

            # Add cursor functionality to the influence line plot
            cursor = mplcursors.cursor(line, hover=True)
            cursor.connect("add", move_arrow)
            cursor.connect("remove", remove_arrow)
            # Update plot
            self.canvas.draw()

    def reset(self):
        self.variable.delete(0, tk.END)
        if hasattr(self, 'spacing_entries'):
            for entry in self.spacing_entries:
                entry.destroy()

        # Clear plots
        self.axes1.clear()
        self.axes2.clear()
        self.axes1.set_title("Beam Diagram")
        self.axes2.set_title("Influence Line")
        self.canvas.draw()

        # Clear cursor annotations
        if self.cursor and hasattr(self.cursor, 'annotations'):
            for annotation in self.cursor.annotations:
                annotation.remove()
            self.cursor.annotations = []

        # Clear the slider and buttons
        if hasattr(self, 'location_slider'):
            self.location_slider.destroy()
        if hasattr(self, 'update_button'):
            self.update_button.destroy()


root = tk.Tk()
print(root.winfo_geometry())

app = BeamInfluenceLine(root)
root.mainloop()