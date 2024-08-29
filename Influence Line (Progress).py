import tkinter as tk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.patches as patches
import numpy as np
import mplcursors
from matplotlib.widgets import Button


class BeamInfluenceLine:
    def __init__(self, root):
        self.root = root
        self.root.title("Influence Line of a Beam")

        # Create input field for number of supports
        self.supports_label = tk.Label(root, text="Number of simple supports (Select between 2 and 5:")
        self.supports_label.pack()
        self.supports_entry = tk.Entry(root)
        self.supports_entry.pack()

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

        self.canvas = FigureCanvasTkAgg(self.figure, master=root)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack()




    def create_spacing_fields(self):
        # Get number of supports

        num_supports = self.supports_entry.get()
        if num_supports:  # Check if num_supports is not empty
            num_supports = int(num_supports)
            # Create input fields for spacing
            self.spacing_entries = []
            #create support names
            self.supports_label=['A','B','C','D','E']

            # Create initial beam diagram for visual purposes
            beam = patches.Rectangle((0, 0.5), (num_supports-1)*5, 2, edgecolor='black', facecolor='blue')
            self.axes1.add_patch(beam)
            self.axes1.arrow(0,4,0,-1,width=.1)
            for i in range(0,(num_supports)):
                if i >= 0:  # Check if spacings list is not empty
                    support = patches.Polygon(np.array(
                        [[i*5, 0.5], [i*5 - .25, -0.5],
                         [i*5 + .25, -0.5]]),
                                              edgecolor='black', facecolor='red')
                    self.axes1.add_patch(support)
                    self.axes1.text(i*5,-1,s=f'{self.supports_label[i]}')
            self.axes1.set_xlim(-1, (num_supports-1)*5 + 1)
            self.axes1.set_ylim(-3, 4)
            #self.axes1.set_aspect('equal', adjustable='box')
            self.canvas.draw()

            for i in range(num_supports-1):
                spacing_label = tk.Label(self.root, text=f"Spacing {self.supports_label[i]} - {self.supports_label[i+1]} (ft):")
                spacing_label.pack()
                spacing_entry = tk.Entry(self.root,width=5,)
                spacing_entry.insert(0,5)

                # Get the position of the subplot in pixels
                #x, y, w, h = self.axes1.get_tk_widget().winfo_x(), self.axes1.get_tk_widget().winfo_y(), self.axes1.get_tk_widget().winfo_width(), self.axes1.get_tk_widget().winfo_height()

                spacing_entry.pack()
                self.spacing_entries.append(spacing_entry)

            # Create button to update plot
            self.update_button = tk.Button(self.root, text="Update Plot", command=self.update_plot)
            self.update_button.pack()



            # Create Slider for influence line section of interest
            self.slider_label = tk.Label(self.root, text="location of influence line ")
            self.slider_label.pack()

            self.location_slider = tk.Scale(self.root, from_=0,
                                            to=1,
                                            resolution=0.01,
                                            orient=tk.HORIZONTAL,
                                            command=self.update_plot, length=900, sliderlength=10, tickinterval=.1)
            self.location_slider.pack()

        else:
            print("Please enter the number of supports")



    def update_plot(self, val=None):
        # Get input values
        num_supports = int(self.supports_entry.get())
        spacings = [float(entry.get()) for entry in self.spacing_entries]
        support_positions = [sum(spacings[:i]) for i in range(num_supports + 1)]
        location = self.location_slider.get()

        # Clear previous plot
        self.axes1.clear()
        self.axes2.clear()
        self.axes1.set_title("Beam Diagram")
        self.axes2.set_title(f'Influence Line @x={round(location*sum(spacings),2)}')
        self.axes1.set_xlabel("Beam Length")
        self.axes2.set_xlabel("Beam Length")
        self.axes2.set_ylabel("Shear influence")

        # Plot beam diagram
        beam = patches.Rectangle((0, 0.5), sum(spacings), 2, edgecolor='black', facecolor='blue')
        self.axes1.add_patch(beam)
        for i, pos in enumerate(support_positions):
            if i <= len(spacings):  # Check if spacings list is not empty
                support = patches.Polygon(np.array([[pos, 0.5], [pos - .25, -0.5], [pos + .25, -0.5]]), edgecolor='black', facecolor='red')
                self.axes1.add_patch(support)
                self.axes1.text(pos, -1, s=f'{self.supports_label[i]}')

        self.axes1.axvline(x=location * sum(spacings), color='r', linestyle='--')  # Add vertical dashed red line
        self.axes1.set_xlim(-1, sum(spacings) + 1)
        self.axes1.set_ylim(-3, 4)
        #self.axes1.set_aspect('equal', adjustable='box')







        # Plot influence line

        if num_supports == 2:
            x = np.arange(0, sum(spacings), .01)
            y = [-i/sum(spacings) if i < location*sum(spacings) else (sum(spacings)-i)/sum(spacings) for i in x]
            self.axes2.text(x=location*sum(spacings),y=-(location*1.2),s= f'{str(round(-location,2))} (max.)')
            self.axes2.text(x=location * sum(spacings), y=(1-location) * 1.2, s=f'{str(round(1-location,2))} (max.)')

        if num_supports == 3:
            x = np.arange(0, sum(spacings), .01)
            # a beam with 3 or more supports is considered statically determinant, thus making the classic equilibrium
            # equations insufficient for any shear calculations. rather, before any influence line equations can be derived,
            # we must first solve a reaction (in terms of the load location) and then can continue
            # with our new adjusted statically "determinant" structure

            # solve the middle support in terms of load location
            L = sum(spacings)  # total length of beam
            c = location*sum(spacings)
            d = L - c
            L1 = spacings[0]  # distance of middle support from LEFT end
            L2 = spacings[1]  # distance of middle support from RIGHT end

            # we will solve the reaction(x) of the middle support by superposition (removing the support
            # and solving the deflection at the middle support location, then solving the required P load at
            # the middle support to produce the same upwards deflection. the following list determines the
            # deflection of the middle support when the load is applied at every position "a" on the beam

            delta_Bx = (L - c) * L1 / (6 * L) * (L ** 2 - (L - c) ** 2 - L1 ** 2) if c >= L1 else c * L2 / (6 * L) * (
                        L ** 2 - c ** 2 - L2 ** 2)
            # this is deflection of middle support as a function of the load location

            # delta_BB= P(L1^2)(L2^2)/3L , this is the deflection @ middle support when upward P load applied at the middle support
            # equating both equations allows us to solve for P (shear @influence location when load placed there):
            B_y = delta_Bx / ((L1 ** 2) * (L2 ** 2) / (3 * L))

            # Now that the middle reaction is determined, the other supports can be solved as functions of load position "a" as well.

            A_y = (-B_y * L2 + (L - c)) / L
            C_y = (-B_y * L1 + c) / L

            #so far we have created values for the reactions of A, B , and C when the load is placed on the beam @ "location"
            #we will use these reaction values to determine the positive and negative influence value at the slider location:
            #using simple statics:

            I_loc_neg=np.array([location,-A_y])
            I_loc_pos=np.array([location,B_y + C_y])

            #now we can add the remaining known Influence values, which is when the load is placed at the supports:
            if location == 0:
                I_Ay=np.array([0,1])
                I_By=np.array([spacings[0],0])
                I_Cy = np.array([spacings[0]+spacings[1],0])
            if location==spacings[0]:
                I_Ay = np.array([0, 0])
                I_By = np.array([spacings[0],1])
                I_Cy = np.array([spacings[0] + spacings[1], 0])
            if location==spacings[0]+spacings[1]:
                I_Ay = np.array([0, 0])
                I_By = np.array([spacings[0],0])
                I_Cy = np.array([spacings[0] + spacings[1], 1])

            else:
                I_Ay = np.array([0, 0])
                I_By = np.array([spacings[0], 0])
                I_Cy = np.array([spacings[0] + spacings[1], 0])





            #we now have 5 datapoints which are known once "location" is given. we can fit these coordinates to a curve
            # by knowing that 2 curves will be needed (before and after "location") and each polynomial will be to the num_supports-1 degree.

            # Extract x and y values
            if location==0 or location==spacings[0] or location==spacings[0]+spacings[1]:
                x1 = np.array([I_Ay[0],I_By[0],I_Cy[0]])
                y1 = np.array([I_Ay[1],I_By[1],I_Cy[1]])
                #create only one polynomial curve if influence location is on support

                # Fit a polynomial curve
                coeffs = np.polyfit(x1, y1, num_supports-1)

                # Generate x values for plotting
                x_plot = np.linspace(min(x1), max(x1), 100)

                # Calculate corresponding y values
                y_plot = np.polyval(coeffs, x_plot)

                x = x_plot
                y = y_plot


            else:
                x1 = np.array([I_Ay[0], I_loc_neg[0], I_loc_pos[0], I_By[0], I_Cy[0]])
                y1 = np.array([I_Ay[1], I_loc_neg[1], I_loc_pos[1], I_By[1], I_Cy[1]])
                print(x1)
                print(y1)


                #create two polynomial curves

                coeffs1 = np.polyfit(x1[:2], y1[:2], num_supports-1)

                # Fit the second polynomial using the remaining data points
                coeffs2 = np.polyfit(x1[2:], y1[2:], num_supports-1)

                # Generate x values for plotting
                x_plot1 = np.linspace(min(x1), x1[1], 100)
                x_plot2 = np.linspace(x1[2], max(x1), 100)

                # Calculate corresponding y values for the first polynomial
                y_plot1 = np.polyval(coeffs1, x_plot1)

                # Calculate corresponding y values for the second polynomial
                y_plot2 = np.polyval(coeffs2, x_plot2)

                x=x_plot1+x_plot2
                y=y_plot1+y_plot2

        #PLACEHOLDER FOR SUPPORTS 4 AND 5
        if num_supports>3:
             x = np.arange(0, sum(spacings), .01)
             y = [1 for i in x]

        # Capture the line object
        line, = self.axes2.plot(x, y)

        self.axes2.plot(x, y)
        self.axes2.axhline(y=0, color='r', linestyle='-')  # Add horizontal red line for x-axis
        self.axes2.set_xlim(-1, sum(spacings)+1)
        self.axes2.set_ylim(-1.5, 1.5)

        # Add cursor functionality to the influence line plot
        cursor = mplcursors.cursor(line, hover=True)
        cursor.connect("add", lambda sel: sel.annotation.set_text(f"x: {sel.target[0]:.2f}, y: {sel.target[1]:.2f}"))

        self.canvas.draw()







root = tk.Tk()
app = BeamInfluenceLine(root)
root.mainloop()
