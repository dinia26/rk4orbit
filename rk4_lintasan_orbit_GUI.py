import numpy as np
import matplotlib.pyplot as plt
from tkinter import *
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

def run_simulation():
    dt = float(dt_entry.get())
    G = float(G_entry.get())
    M = float(M_entry.get())
    npoints = int(npoints_entry.get())

    x0 = (float(planet_distances[planet_var.get()])) * 1.496e11
    y = 0
    R = np.sqrt(x0**2 + y**2)

    R_s = (float(planet_average_distance[planet_var.get()])) * 1.496e11

    s = ((float(planet_a[planet_var.get()]) + float(planet_b[planet_var.get()])) / 2) * 1.496e11
    
    mass = planet_masses[planet_var.get()]
    
    v0 = np.sqrt((G * M / s) * ((2 * s / R) - 1))       #ubah variasi kecepatan sesuai dengan yang diinginkan
    v0_entry.delete(0, END)
    v0_entry.insert(0, "{:.6f}".format(v0))
    
    vx = 0
    vy = v0
    
    x, y = x0, 0
    x_values = [x]
    y_values = [y]
    r_values = [R]

    vx_half_period = 0
    vy_half_period = 0
    
    orbit_steps = None
    orbit_time = None

    aphelion = x0
    perihelion = x0

    for step in range(1, npoints):                      #metode: Runge Kutta Orde 4
        k1x = dt * vx
        k1y = dt * vy

        R = np.sqrt((x + k1x/2)**2 + (y + k1y/2)**2)

        k1vx = dt * (-G * M / R**3) * x
        k1vy = dt * (-G * M / R**3) * y

        k2x = dt * (vx + k1vx / 2)
        k2y = dt * (vy + k1vy / 2)

        k2vx = dt * (-G * M / R**3 * (x + k1x / 2))
        k2vy = dt * (-G * M / R**3 * (y + k1y / 2))

        k3x = dt * (vx + k2vx / 2)
        k3y = dt * (vy + k2vy / 2)

        k3vx = dt * (-G * M / R**3 * (x + k2x / 2))
        k3vy = dt * (-G * M / R**3 * (y + k2y / 2))

        k4x = dt * (vx + k3vx)
        k4y = dt * (vy + k3vy)

        k4vx = dt * (-G * M / R**3 * (x + k3x))
        k4vy = dt * (-G * M / R**3 * (y + k3y))    

        x = x + 1/6 * (k1x + 2 * k2x + 2 * k3x + k4x)
        y = y + 1/6 * (k1y + 2 * k2y + 2 * k3y + k4y)

        vx = vx + 1/6 * (k1vx + 2 * k2vx + 2 * k3vx + k4vx)
        vy = vy + 1/6 * (k1vy + 2 * k2vy + 2 * k3vy + k4vy)

        x_values.append(x)
        y_values.append(y)
        r_values.append(R)

        if y < 0 and y_values[step - 1] >= 0 and orbit_steps is None:
            orbit_steps = step
            orbit_time = step * dt

    if orbit_steps is not None:
        vx_half_period = vx
        vy_half_period = vy
        #half_period_velocity_entry.delete(0, END)
        #half_period_velocity_entry.insert(0, "vx: {:.6f}, vy: {:.6f}".format(vx_half_period, vy_half_period))

    if R > aphelion:
        aphelion = R
    if R < perihelion:
        perihelion = R

    if orbit_time is not None:
        period = 2 * orbit_time / (24*3600)
        period_entry.delete(0, END)
        period_entry.insert(0, "{:.6f}".format(period))
    else:
        period_entry.delete(0, END)
        period_entry.insert(0, "Tidak dapat menghitung")

    semi_major_axis = ((max(x_values) + abs(min(x_values))) / 2) / 1e3
    semi_minor_axis = ((max(y_values) + abs(min(y_values))) / 2) / 1e3
    semi_major_entry.delete(0, END)
    semi_major_entry.insert(0, "{:.6f}".format(semi_major_axis))
    semi_minor_entry.delete(0, END)
    semi_minor_entry.insert(0, "{:.6f}".format(semi_minor_axis))

    aphelion = max(r_values) / 1e3
    perihelion = min(r_values) / 1e3
    eccentricity = (aphelion - perihelion) / (aphelion + perihelion)
    eccentricity_entry.delete(0, END)
    eccentricity_entry.insert(0, "{:.13f}".format(eccentricity))

    R_val = ((aphelion + perihelion) / 2) * 1e3
    
    if orbit_steps is not None:
        PE = -G * M * mass / R_val

        v_half = np.sqrt(vx_half_period**2 + vy_half_period**2)

        v_R = (abs(v_half) + abs(v0)) / 2
        
        #KE_x = 0.5 * mass * vx**2
        #KE_y = 0.5 * mass * vy**2

        KE_total = 0.5 * mass * v_R**2

        E_total = KE_total + PE
        potential_energy_entry.delete(0, END)
        potential_energy_entry.insert(0, "{:.10f}".format(PE))
        kinetic_energy_entry.delete(0, END)
        kinetic_energy_entry.insert(0, "{:.10f}".format(KE_total))
        total_energy_entry.delete(0, END)
        total_energy_entry.insert(0, "{:.10f}".format(E_total))
    else:
        potential_energy_entry.delete(0, END)
        potential_energy_entry.insert(0, "Tidak dapat menghitung")
        kinetic_energy_entry.delete(0, END)
        kinetic_energy_entry.insert(0, "Tidak dapat menghitung")
        total_energy_entry.delete(0, END)
        total_energy_entry.insert(0, "Tidak dapat menghitung")

    ax.clear()
    ax.plot(0, 0, 'oy', markersize=8, markerfacecolor='yellow')
    ax.plot(x0, 0, 'ro')
    ax.plot(x_values, y_values)
    if orbit_steps is not None:
        ax.plot(x_values[orbit_steps - 1], y_values[orbit_steps - 1], 'ro', label='Periode (1/2 putaran)')
        ax.text(x_values[orbit_steps - 1], y_values[orbit_steps - 1], 'Periode (1/2 putaran)', fontsize=10, color='red', verticalalignment='bottom')
    ax.axis('equal')
    ax.set_xlabel('x (m)', fontsize=14)  
    ax.set_ylabel('y (m)', fontsize=14)  
    ax.grid(True)
    ax.set_title('Simulasi Orbit Planet {}'.format(planet_var.get()), fontsize=16)  
    canvas.draw()

    perihelion_entry.delete(0, END)
    perihelion_entry.insert(0, "{:.10f}".format(perihelion))
    aphelion_entry.delete(0, END)
    aphelion_entry.insert(0, "{:.10f}".format(aphelion))

# GUI setup
root = Tk()
root.title("Simulasi Orbit Planet")

input_frame = Frame(root)
input_frame.pack(side=LEFT, padx=20, pady=20)

plot_frame = Frame(root)
plot_frame.pack(side=LEFT, padx=20, pady=20) 

result_frame = Frame(root)
result_frame.pack(side=LEFT, padx=20, pady=20)

planet_var = StringVar()
planet_var.set("Merkurius")
Label(input_frame, text="Pilih Planet:", font=("Helvetica", 14)).pack()   
planet_menu = OptionMenu(input_frame, planet_var, "Merkurius", "Venus", "Bumi", "Mars", "Jupiter", "Saturnus", "Uranus", "Neptunus")
planet_menu.config(font=("Helvetica", 14))   
planet_menu.pack()

planet_distances = {
    "Merkurius": 0.31,
    "Venus": 0.72,
    "Bumi": 0.98,
    "Mars": 1.38,
    "Jupiter": 4.95,
    "Saturnus": 9.01,
    "Uranus": 18.28,
    "Neptunus": 29.8
}

planet_masses = {
    "Merkurius": 0.330e24,
    "Venus": 4.87e24,
    "Bumi": 5.97e24,
    "Mars": 0.642e24,
    "Jupiter": 1898e24,
    "Saturnus": 568e24,
    "Uranus": 86.8e24,
    "Neptunus": 102e24
}

planet_a = {
    "Merkurius": 0.47,
    "Venus": 0.73,
    "Bumi": 1.02,
    "Mars": 1.67,
    "Jupiter": 5.45,
    "Saturnus": 10.07,
    "Uranus": 20.09,
    "Neptunus": 30.32
}

planet_b = {
    "Merkurius": 0.31,
    "Venus": 0.72,
    "Bumi": 0.98,
    "Mars": 1.38,
    "Jupiter": 4.95,
    "Saturnus": 9.01,
    "Uranus": 18.28,
    "Neptunus": 29.8
}

planet_average_distance = {
    "Merkurius": 0.387,
    "Venus": 0.723,
    "Bumi": 1,
    "Mars": 1.523,
    "Jupiter": 5.202,
    "Saturnus": 9.539,
    "Uranus": 19.18,
    "Neptunus": 30.06
}

Label(input_frame, text="dt:", font=("Helvetica", 14)).pack()   
dt_entry = Entry(input_frame)
dt_entry.config(font=("Helvetica", 14))   
dt_entry.pack()
dt_entry.insert(0, "10000")

Label(input_frame, text="G:", font=("Helvetica", 14)).pack()   
G_entry = Entry(input_frame)
G_entry.config(font=("Helvetica", 14))   
G_entry.pack()
G_entry.insert(0, "0.000000000066743")

Label(input_frame, text="M:", font=("Helvetica", 14)).pack()   
M_entry = Entry(input_frame)
M_entry.config(font=("Helvetica", 14))   
M_entry.pack()
M_entry.insert(0, "1988400000000000000000000000000")

Label(input_frame, text="npoints:", font=("Helvetica", 14)).pack()   
npoints_entry = Entry(input_frame)
npoints_entry.config(font=("Helvetica", 14))   
npoints_entry.pack()
npoints_entry.insert(0, "1000000")

Label(input_frame, text="v0 (m/s):", font=("Helvetica", 14)).pack()   
v0_entry = Entry(input_frame)
v0_entry.config(font=("Helvetica", 14))
v0_entry.pack()

Label(input_frame, text="Perihelion:", font=("Helvetica", 14)).pack()   
perihelion_entry = Entry(input_frame)
perihelion_entry.config(font=("Helvetica", 14))
perihelion_entry.pack()

Label(input_frame, text="Aphelion:", font=("Helvetica", 14)).pack()   
aphelion_entry = Entry(input_frame)
aphelion_entry.config(font=("Helvetica", 14))
aphelion_entry.pack()

Button(input_frame, text="Run Simulation", command=run_simulation, font=("Helvetica", 14)).pack()

fig, ax = plt.subplots(figsize=(8, 6))  
canvas = FigureCanvasTkAgg(fig, master=plot_frame)
canvas.get_tk_widget().pack()
canvas.draw()

Label(result_frame, text="Periode orbit: ", font=("Helvetica", 14)).pack()
period_entry = Entry(result_frame, font=("Helvetica", 14))
period_entry.pack()

Label(result_frame, text="Sumbu Semi Mayor: ", font=("Helvetica", 14)).pack()
semi_major_entry = Entry(result_frame, font=("Helvetica", 14))
semi_major_entry.pack()

Label(result_frame, text="Sumbu Semi Minor: ", font=("Helvetica", 14)).pack()
semi_minor_entry = Entry(result_frame, font=("Helvetica", 14))
semi_minor_entry.pack()

Label(result_frame, text="Eksentrisitas orbit: ", font=("Helvetica", 14)).pack()
eccentricity_entry = Entry(result_frame, font=("Helvetica", 14))
eccentricity_entry.pack()

Label(result_frame, text="Energi Potensial: ", font=("Helvetica", 14)).pack()
potential_energy_entry = Entry(result_frame, font=("Helvetica", 14))
potential_energy_entry.pack()

Label(result_frame, text="Energi Kinetik Total: ", font=("Helvetica", 14)).pack()
kinetic_energy_entry = Entry(result_frame, font=("Helvetica", 14))
kinetic_energy_entry.pack()

Label(result_frame, text="Energi Total: ", font=("Helvetica", 14)).pack()
total_energy_entry = Entry(result_frame, font=("Helvetica", 14))
total_energy_entry.pack()

#Label(result_frame, text="Kecepatan di Setengah Periode: ", font=("Helvetica", 14)).pack()
#half_period_velocity_entry = Entry(result_frame, font=("Helvetica", 14))
#half_period_velocity_entry.pack()

root.mainloop()
