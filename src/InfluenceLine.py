import tkinter as tk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.patches as patches
import numpy as np
import mplcursors


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


        #PLACEHOLDER FOR SUPPORTS 4 AND 5
        if num_supports>2:
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


        # Update plot
        self.canvas.draw()



root = tk.Tk()
app = BeamInfluenceLine(root)
root.mainloop()

