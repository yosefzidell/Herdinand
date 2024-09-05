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
        self.root.title("Influence Line")

        # Reset Button
        self.reset_button = tk.Button(root, text="Reset", command=self.reset)
        self.reset_button.place(x=.9,y=.5)

        # Create input field for number of supports
        self.supports_label = tk.Label(root, text="step 1: click on dropdown menu below to select number of supports.")
        self.supports_label.pack()
        #
        def create_dropdown(options, master):
            variable = tk.StringVar(master)
            variable.set(options[0])  # default value
            dropdown = tk.OptionMenu(master, variable, *options)
            return dropdown, variable

        self.options = [2, 3, 4]
        self.dropdown, self.variable = create_dropdown(self.options, root)

        self.dropdown.pack()

        self.force_label = tk.Label(root, text="step 2: click on dropdown menu below to select type of influence line desired. then press next.")
        self.force_label.pack()

        def create_dropdown2(forces, master):
            variable2 = tk.StringVar(master)
            variable2.set(forces[0])  # default value
            dropdown2 = tk.OptionMenu(master, variable2, *forces)
            return dropdown2, variable2

        self.forces = ['Shear Influence','Moment Influence']
        self.dropdown2, self.variable2 = create_dropdown2(self.forces, root)

        self.dropdown2.pack()


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
        self.bubble_S1 = (plt.Rectangle((0.09, .15), 0.01,height=.05, facecolor='blue', edgecolor='black', transform=self.axes1.transAxes, visible=False))
        self.axes1.add_patch(self.bubble_S1)
        self.bubble_S1_status = True
        self.bubble_M1 = (plt.Rectangle((0.09, .05), 0.01,height=.05, facecolor='white', edgecolor='black', transform=self.axes1.transAxes, visible=False))
        self.axes1.add_patch(self.bubble_M1)
        self.bubble_M1_status = False
        self.bubble_S2 = (plt.Rectangle((.92, .15), 0.01, height=.05, facecolor='blue', edgecolor='black',
                                       transform=self.axes1.transAxes, visible=False))
        self.bubble_S2_status = True
        self.axes1.add_patch(self.bubble_S2)
        self.bubble_M2 = (plt.Rectangle((.92, .05), 0.01, height=.05, facecolor='white', edgecolor='black',
                                        transform=self.axes1.transAxes, visible=False))
        self.axes1.add_patch(self.bubble_M2)
        self.bubble_M2_status = False

        self.canvas = FigureCanvasTkAgg(self.figure, master=root)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack()

        self.supports_options=[self.bubble_S1,self.bubble_M1,self.bubble_S2,self.bubble_M2]


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
        self.info = tk.Label(root, text="step 3: determine the support types in the beam diagram, and the spacing between supports below. then click 'update plot'. ")
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
            support = patches.Polygon(np.array(
                [[i * 5, 0.5], [i * 5 - .25, -0.5],
                 [i * 5 + .25, -0.5]]),
                edgecolor='black', facecolor='red')
            self.axes1.add_patch(support)
            self.axes1.text(i * 5, -1, s=f'{self.supports_label[i]}')
        self.axes1.set_xlim(-1, (num_supports - 1) * 5 + 1)
        self.axes1.set_ylim(-3, 4)


        # make selection bubbles visible, move them to be directly under the end supports, and add in text
        self.axes1.text(.07, .15, 'S',color='blue', transform=self.axes1.transAxes)
        self.bubble_S1.set_visible(True)
        self.axes1.text(.07, .05, 'M',color='blue', transform=self.axes1.transAxes)
        self.bubble_M1.set_visible(True)
        self.axes1.text(.9, .15, 'S', color='blue', transform=self.axes1.transAxes)
        self.bubble_S2.set_visible(True)
        self.axes1.text(.9, .05, 'M', color='blue', transform=self.axes1.transAxes)
        self.bubble_M2.set_visible(True)

        self.canvas.draw()
        self.canvas.mpl_connect('button_press_event', self.on_pick)

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
        self.slider_label = tk.Label(self.root, text="step 4: drag slider below to the influence position on the beam, indicated by the dashed red line in the Beam Diagram. you can hover over the influence Line to see the influence value.")
        self.slider_label.pack()

        self.location_slider = tk.Scale(self.root, from_=0,
                                        to=1,
                                        resolution=0.01,
                                        orient=tk.HORIZONTAL,
                                        command=self.update_plot, length=900, sliderlength=10, tickinterval=.1)
        self.location_slider.pack()




    # Function to handle support click

    def on_pick(self,event):
        if event.inaxes == self.axes1:
            if self.bubble_S1.contains(event)[0]:
                self.bubble_S1.set_facecolor('blue')
                self.bubble_M1.set_facecolor('white')
                self.bubble_S1_status = True
                self.bubble_M1_status = False
            elif self.bubble_M1.contains(event)[0]:
                self.bubble_S1.set_facecolor('white')
                self.bubble_M1.set_facecolor('blue')
                self.bubble_S1_status = False
                self.bubble_M1_status = True
            elif self.bubble_S2.contains(event)[0]:
                self.bubble_S2.set_facecolor('blue')
                self.bubble_M2.set_facecolor('white')
                self.bubble_S2_status = True
                self.bubble_M2_status = False
            elif self.bubble_M2.contains(event)[0]:
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
            self.axes2.set_title(f'{self.variable2.get()} line @x={round(location * sum(spacings), 2)}')
            self.axes1.set_xlabel("Beam Length (ft)")
            self.axes2.set_xlabel("Beam Length (ft)")
            self.axes2.set_ylabel(f'{self.variable2.get()}')

            # Plot beam diagram
            beam = patches.Rectangle((0, 0.5), sum(spacings), 2, edgecolor='black', facecolor='grey')
            self.axes1.add_patch(beam)
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
                if self.bubble_S1_status==True and self.bubble_S2_status==True:
                    x = np.arange(0, sum(spacings), .01)
                    L = sum(spacings)
                    c = location * L
                    A_y=[i/L for i in x]
                    if self.variable2.get()==self.forces[0]: #shear influence line selected
                     y = [-i if j < c else -i+1 for i,j in zip(A_y,x)]
                    else: #moment influence line selected
                     y = [(i)*(L-c)/L if i<=c else (L-c)*c/L-(i-c)*((L-c)*c/L)/(L-c) for i in x]

                elif self.bubble_M1_status == True and self.bubble_S2_status == True:
                    #one end is a moment connection and the other is a simple support:
                    x = np.arange(0, sum(spacings), .01)
                    #remove support B and solve its deflection
                    L=sum(spacings)
                    c=location*L
                    delta_Bx=np.array([i**2/6*(3*L-i) for i in x])
                    #solve reaction By with adding the support back in thats equal to delta_Bx
                    B_y = np.array([i*3/(L**3) for i in delta_Bx])
                    #now that B_y is solved, solve M_A and A_y usng equilibrium:
                    M_A=np.array([-i*L+j for i,j in zip(B_y,x)])
                    A_y=np.array([(i+L-j)/L for i,j in zip(M_A,x)])

                    # now we can determine the shear influence using A_y
                    y = np.array([i - 1 if j < c else i for i, j in zip(A_y, x)])

                elif self.bubble_S1_status == True and self.bubble_M2_status == True:
                    #one end is a moment connection and the other is a simple support:
                    x = np.arange(0, sum(spacings), .01)
                    #remove support A and solve its deflection
                    L=sum(spacings)
                    c=location*L
                    delta_Ax=np.array([(L-i)**2/6*(3*L-(L-i)) for i in x])
                    #solve reaction By with adding the support back in thats equal to delta_Bx
                    A_y = np.array([i*3/(L**3) for i in delta_Ax])
                    #now that B_y is solved, solve M_A and A_y usng equilibrium:
                    M_B=np.array([-i*L+(L-j) for i,j in zip(A_y,x)])
                    B_y=np.array([(i+j)/L for i,j in zip(M_B,x)])

                    y = np.array([i - 1 if j < c else i for i, j in zip(A_y, x)])





                else:
                    # both ends are moments supports
                    x = np.arange(0, sum(spacings), .01)
                    # solve for right side support Ay and M_A
                    L = sum(spacings)
                    c = location * L

                    A_y = np.array([((L-i)**2/(L**3))*(3*i+(L-i)) for i in x])
                    M_A = np.array([(i*(L-i)**2)/(L**2) for i in x])

                    # now we can determine the shear influence using A_y and M_A
                    # now we can determine the influence using A_y and M_A
                    if self.variable2.get() == self.forces[0]:  # shear influence line selected
                        y = [i-1 if j <= c else i for i, j in zip(A_y, x)]
                    else:  # moment influence line selected
                        y = [i * c - k -(c-j) if j <= c else i * c - k for i, j, k in zip(A_y, x, M_A)]

            if num_supports == 3:
                x = np.arange(0, sum(spacings), .01)
                L = sum(spacings)  # total length of beam
                c = L * location  # location of interest for influence line diagram, determined by the slider
                L1 = spacings[0]  # distance of middle support from left end
                L2 = spacings[1]  # distance of middle support from right end

                if self.bubble_S1_status == True and self.bubble_S2_status == True:
                    # a beam with 3 or more supports is considered statically indeterminant, thus making the classic equilibrium
                    # equations insufficient for any shear calculations. rather, before any influence line equations can be derived,
                    # we must first solve a reaction in terms of the load location "x" and then can continue
                    # with our new adjusted statically "determinant" structure.



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
                    if c < L1:
                        y = np.array([i - 1 if j < c else i for i, j in zip(A_y, x)])
                    elif c == L1:
                        y = np.array([i for i in B_y])
                    elif c<L:
                        y = np.array([(i + k - 1) if j < c else (i + k) for i, k, j in zip(A_y, B_y, x)])
                    else:
                        y = np.array([-i for i in C_y])
                elif self.bubble_S1_status == True and self.bubble_M2_status == True:
                    #Moment support on right
                    # determine middle support By using its redundant deflection
                    delta_Bx = np.array(
                        [(L - i)**2 * L1 / (12 * L**3) * (3*i*L**2 - 2*L*L1**2 - i*L1 ** 2) if i >= L1 else i / (12 * L**3) * (L2**2)*(
                                3*L ** 2*L1 - L1*i ** 2 - 2*L*i ** 2) for i in x])
                    B_y=np.array([(i*12*L**3)/((L1**2*L2**3)*(3*L+L1)) for i in delta_Bx])

                    #solve left support A_y by using its redundant deflection. remove support A, solve deflection at
                    #free end from Unit Load and B_y, and solve A_y by making it equal that deflection.

                    delta_B_cantilever=np.array([L2**2/6*(3*(L-i)-L+L1) if i<=L1 else (L-i)**2/6*(3*L-3*L1-(L-i)) for i in x])
                    #B_y_cantilever=np.array([3*i/((L-j)**3) for i,j in zip(delta_B_cantilever,x)])
                    delta_A1=np.array([(L-i)**2/6*(3*L-(L-i)) for i in x])
                    delta_A2=np.array([i*(L2**2)/6*(3*L-L2) for i in B_y])
                    A_y=np.array([3*(i-j)/(L**3) for i,j in zip(delta_A1,delta_A2)])
                    #now y can be determined by using A_y and B_y
                    if self.variable2.get() == self.forces[0]:  # shear influence line selected
                        if c < L1:
                            y = np.array([i - 1 if j < c else i for i, j in zip(A_y, x)])
                        elif c == L1:
                            y = np.array([i for i in B_y])
                        else:
                            y = np.array([(i + k - 1) if j < c else (i + k) for i, k, j in zip(A_y, B_y, x)])
                    else:  # moment influence line selected
                        if c <= L1:
                            y = np.array([i*c -(c-j) if j < c else i*c for i, j in zip(A_y, x)])
                        else:
                            y = np.array([i*c -(c-j)+k*(c-L1) if j < c else i*c +k*(c-L1) for i, j, k in zip(A_y, x,B_y)])
                elif self.bubble_M1_status == True and self.bubble_S2_status == True:
                    #Moment support on left
                    #the influence line will mirror the previous condition , thus the same procedure will be used by changing x to (L-x) and L1 to L2, and vice versa
                    # determine middle support By using its redundant deflection
                    delta_Bx = np.array(
                        [(i)**2 * L2 / (12 * L**3) * (3*(L-i)*L**2 - 2*L*L2**2 - (L-i)*L1 ** 2) if (L-i) >= L2 else (L-i) / (12 * L**3) * (L1**2)*(
                                3*L ** 2*L2 - L2*(L-i) ** 2 - 2*L*(L-i) ** 2) for i in x])
                    B_y=np.array([(i*12*L**3)/((L2**2*L1**3)*(3*L+L2)) for i in delta_Bx])

                    #solve right support A_y by using its redundant deflection. remove support A, solve deflection at
                    #free end from Unit Load and B_y, and solve A_y by making it equal that deflection.

                    delta_B_cantilever=np.array([L2**2/6*(3*(L-i)-L+L1) if i<=L1 else (L-i)**2/6*(3*L-3*L1-(L-i)) for i in x])
                    #B_y_cantilever=np.array([3*i/((L-j)**3) for i,j in zip(delta_B_cantilever,x)])
                    delta_A1=np.array([(i)**2/6*(3*L-(i)) for i in x])
                    delta_A2=np.array([i*(L1**2)/6*(3*L-L1) for i in B_y])
                    A_y=np.array([3*(i-j)/(L**3) for i,j in zip(delta_A1,delta_A2)])
                    #now y can be determined by using A_y and B_y
                    if self.variable2.get() == self.forces[0]:  # shear influence line selected
                        if (L-c) < L2:
                            y = np.array([i - 1 if (L-j) < (L-c) else i for i, j in zip(A_y, x)])
                        elif (L-c) == L2:
                            y = np.array([i for i in B_y])
                        else:
                            y = np.array([(i + k - 1) if (L-j) < (L-c) else (i + k) for i, k, j in zip(A_y, B_y, x)])
                    else:  # moment influence line selected
                        if (L-c) <= L2:
                            y = np.array([i*(L-c) -((L-c)-(L-j)) if (L-j) < (L-c) else i*(L-c) for i, j in zip(A_y, x)])
                        else:
                            y = np.array([i*(L-c) -((L-c)-(L-j))+k*((L-c)-L2) if (L-j) < (L-c) else i*(L-c) +k*((L-c)-L2) for i, j, k in zip(A_y, x,B_y)])
                else:
                    # both ends are moments supports
                    x = np.arange(0, sum(spacings), .01)
                    delta_Bx = np.array(
                        [(L - i)**2*L1**2 / ( 6 * L ** 3) * (3 * i*L - 3*L1 * i - L1 * (L - i)) if i>L1 else (i) ** 2 * L2**2 / (6 * L ** 3) * (3 * (L - i) * L - 3 * (L-i) * L2 - i * L2 ) for i in x])
                    B_y = np.array([(i * 3 * L**3) / (L2 ** 3 * L1 ** 3) for i in delta_Bx])

                    #now solve for M_A and A_y by determining M_A an A_y from the effects of the unit load, and B_y individually.

                    A_y1 = np.array([((L-i)**2/(L**3))*(3*i+(L-i)) for i in x]) #A_y from unit load
                    M_A1 = np.array([(i*(L-i)**2)/(L**2) for i in x]) #M_A from unit load

                    A_y2 = np.array([(i*(L2) ** 2 / (L ** 3)) * (3 * L1 + L2) for i in B_y])  # A_y from B_y
                    M_A2 = np.array([(i * (L1)*(L2) ** 2) / (L ** 2) for i in B_y])  # M_A from B_y

                    A_y=np.array([i-j for i,j in zip(A_y1,A_y2)])
                    M_A=np.array([i-j for i,j in zip(M_A2,M_A1)])


                    # now we can determine the influence using A_y and M_A
                    if self.variable2.get() == self.forces[0]:  # shear influence line selected
                        if c==0:
                            y = np.array([i for i in A_y])
                        if c < L1:
                            y = np.array([i - 1 if j < c else i for i, j in zip(A_y, x)])
                        elif c == L1:
                            y = np.array([i for i in B_y])
                        else:
                            y = np.array([(i + k - 1) if j < c else (i + k) for i, k, j in zip(A_y, B_y, x)])
                    else:  # moment influence line selected
                        if c <= L1:
                            y = np.array([i * c + k - (c - j) if j < c else i * c +k for i, j,k in zip(A_y, x,M_A)])
                        else:
                            y = np.array(
                                [i * c + k - (c - j)+b*(c-L1) if j < c else i * c +k +b*(c-L1) for i, j,k,b in zip(A_y, x,M_A,B_y)])







            # PLACEHOLDER FOR SUPPORTS 4 AND 5
            if num_supports > 3:
                x = np.arange(0, sum(spacings), .01)
                L = sum(spacings)  # total length of beam
                c = L * location  # location of interest for influence line diagram, determined by the slider
                y = [1 for i in x]
                # Add text to the subplot
                self.axes2.text(L/3,1.25, 'Influence Line not available yet', fontsize=15, color='red')


            # Capture the line object
            line = self.axes2.plot(x, y)

            self.axes2.plot(x, y)
            self.axes2.axhline(y=0, color='r', linestyle='--')  # Add horizontal red line for x-axis
            self.axes2.set_xlim(-1, sum(spacings) + 1)
            self.axes2.set_ylim(min(y)-1, max(y)+1)

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