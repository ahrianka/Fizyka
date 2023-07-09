import scipy as sci
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
import matplotlib
matplotlib.use('TkAgg')  # Use the TkAgg backend


G = 6.67408e-11  # N-m2/kg2

m_nd = 1.989e+30  # kg #mass of the sun
r_nd = 5.326e+12  # m #distance between stars in Alpha Centauri
v_nd = 30000  # m/s #relative velocity of earth around the sun
t_nd = 79.91 * 365 * 24 * 3600 * 0.51  # s #orbital period of Alpha Centauri

K1 = G * t_nd * m_nd / (r_nd ** 2 * v_nd)
K2 = v_nd * t_nd / r_nd


# masses
m1 = 1.1  # Alpha Centauri A
m2 = 0.9  # Alpha Centauri B
m3 = 1.0  # Third Star

# m1 = 0.5  # Alpha Centauri A
# m2 = 1.5  # Alpha Centauri B
# m3 = 0.1  # Third Star
#
# m1 = 1.0  # Alpha Centauri A
# m2 = 0.01  # Alpha Centauri B
# m3 = 0.5  # Third Star


# position vectors
r1 = [-0.5, 0, 0]  # m
r2 = [0.5, 0, 0]  # m
r3 = [0, 1, 0]  # m

# r1 = [-0.5, 0, 0]  # m
# r2 = [0.5, 0, 0]  # m
# r3 = [-0.5, 0, -0.5]  # m

r1 = sci.array(r1, dtype="float64")
r2 = sci.array(r2, dtype="float64")
r3 = sci.array(r3, dtype="float64")

# centr masy
r_com = (m1 * r1 + m2 * r2) / (m1 + m2)

# speed m/s
v1 = [0.01, 0.01, 0]
v2 = [-0.05, 0, -0.1]
v3 = [0, -0.01, 0]

# to array
v1 = sci.array(v1, dtype="float64")
v2 = sci.array(v2, dtype="float64")
v3 = sci.array(v3, dtype="float64")

# prędkość centru masy
v_com = (m1 * v1 + m2 * v2) / (m1 + m2)



def TwoBodyEquations(w, t, G, m1, m2):
    r1 = w[:3]
    r2 = w[3:6]
    v1 = w[6:9]
    v2 = w[9:12]
    r = sci.linalg.norm(r2 - r1)
    dv1bydt = K1 * m2 * (r2 - r1) / r ** 3
    dv2bydt = K1 * m1 * (r1 - r2) / r ** 3
    dr1bydt = K2 * v1
    dr2bydt = K2 * v2
    r_derivs = sci.concatenate((dr1bydt, dr2bydt))
    derivs = sci.concatenate((r_derivs, dv1bydt, dv2bydt))
    return derivs



init_params = sci.array([r1, r2, v1, v2])  # create array of initial params
init_params = init_params.flatten()  # flatten array to make it 1D
time_span = sci.linspace(0, 8, 500)  # 8 orbital periods and 500 points


two_body_sol = sci.integrate.odeint(TwoBodyEquations, init_params, time_span, args=(G, m1, m2))

r1_sol = two_body_sol[:, :3]
r2_sol = two_body_sol[:, 3:6]


fig = plt.figure(figsize=(15, 15))

ax = fig.add_subplot(111, projection="3d")

ax.plot(r1_sol[:, 0], r1_sol[:, 1], r1_sol[:, 2], color="darkblue")
ax.plot(r2_sol[:, 0], r2_sol[:, 1], r2_sol[:, 2], color="tab:red")

ax.scatter(r1_sol[-1, 0], r1_sol[-1, 1], r1_sol[-1, 2], color="darkblue", marker="o", s=100, label="Alpha Centauri A")
ax.scatter(r2_sol[-1, 0], r2_sol[-1, 1], r2_sol[-1, 2], color="tab:red", marker="o", s=100, label="Alpha Centauri B")

ax.set_xlabel("x-coordinate", fontsize=14)
ax.set_ylabel("y-coordinate", fontsize=14)
ax.set_zlabel("z-coordinate", fontsize=14)
ax.set_title("Visualization of orbits of stars in a two-body system\n", fontsize=14)
ax.legend(loc="upper left", fontsize=14)

rcom_sol = (m1 * r1_sol + m2 * r2_sol) / (m1 + m2)
r1com_sol = r1_sol - rcom_sol
r2com_sol = r2_sol - rcom_sol



r_com = (m1 * r1 + m2 * r2 + m3 * r3) / (m1 + m2 + m3)
v_com = (m1 * v1 + m2 * v2 + m3 * v3) / (m1 + m2 + m3)


def ThreeBodyEquations(w, t, G, m1, m2, m3):
    r1 = w[:3]
    r2 = w[3:6]
    r3 = w[6:9]
    v1 = w[9:12]
    v2 = w[12:15]
    v3 = w[15:18]
    r12 = sci.linalg.norm(r2 - r1)
    r13 = sci.linalg.norm(r3 - r1)
    r23 = sci.linalg.norm(r3 - r2)

    dv1bydt = K1 * m2 * (r2 - r1) / r12 ** 3 + K1 * m3 * (r3 - r1) / r13 ** 3
    dv2bydt = K1 * m1 * (r1 - r2) / r12 ** 3 + K1 * m3 * (r3 - r2) / r23 ** 3
    dv3bydt = K1 * m1 * (r1 - r3) / r13 ** 3 + K1 * m2 * (r2 - r3) / r23 ** 3
    dr1bydt = K2 * v1
    dr2bydt = K2 * v2
    dr3bydt = K2 * v3
    r12_derivs = sci.concatenate((dr1bydt, dr2bydt))
    r_derivs = sci.concatenate((r12_derivs, dr3bydt))
    v12_derivs = sci.concatenate((dv1bydt, dv2bydt))
    v_derivs = sci.concatenate((v12_derivs, dv3bydt))
    derivs = sci.concatenate((r_derivs, v_derivs))
    return derivs


init_params = sci.array([r1, r2, r3, v1, v2, v3])
init_params = init_params.flatten()
time_span = sci.linspace(0, 20, 500)

three_body_sol = sci.integrate.odeint(ThreeBodyEquations, init_params, time_span, args=(G, m1, m2, m3))

r1_sol = three_body_sol[:, :3]
r2_sol = three_body_sol[:, 3:6]
r3_sol = three_body_sol[:, 6:9]

fig = plt.figure(figsize=(15, 15))
ax = fig.add_subplot(111, projection="3d")


#animacja
def update_animation(i):
    ax.clear()
    # Plot the orbits
    ax.plot(r1_sol[:i, 0], r1_sol[:i, 1], r1_sol[:i, 2], color="darkblue")
    ax.plot(r2_sol[:i, 0], r2_sol[:i, 1], r2_sol[:i, 2], color="tab:red")
    ax.plot(r3_sol[:i, 0], r3_sol[:i, 1], r3_sol[:i, 2], color="green")
    # Plot the final positions of the stars
    ax.scatter(r1_sol[i, 0], r1_sol[i, 1], r1_sol[i, 2], color="darkblue", marker="o", s=100, label="Alpha Centauri A")
    ax.scatter(r2_sol[i, 0], r2_sol[i, 1], r2_sol[i, 2], color="tab:red", marker="o", s=100, label="Alpha Centauri B")
    ax.scatter(r3_sol[i, 0], r3_sol[i, 1], r3_sol[i, 2], color="green", marker="o", s=100, label="Third Star")
    # Add labels and title
    ax.set_xlabel("x-coordinate", fontsize=14)
    ax.set_ylabel("y-coordinate", fontsize=14)
    ax.set_zlabel("z-coordinate", fontsize=14)
    ax.set_title("Visualization of orbits of stars in a three-body system\n", fontsize=14)
    ax.legend(loc="upper left", fontsize=14)



ani = animation.FuncAnimation(fig, update_animation, frames=len(time_span), interval=50)
fig.set_size_inches(10, 6)

plt.show()
