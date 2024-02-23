import sys
# Add paths
sys.path.insert(1, './GNC')
sys.path.insert(1, './graphics')
sys.path.insert(1, './postproc')
sys.path.insert(1, './vehicles')
sys.path.insert(1, './algorithm')

import agents as ag
import pygame
import matplotlib.pyplot as pl
import matplotlib.colors as mcolors
import drawmisc
import agents as ag
import gradiente as grad
import numpy as np
from numpy import linalg as la

import logpostpro as lp
import gvf

# setup simulation
WIDTH = 1200
HEIGHT = 700

CENTERX = WIDTH/2
CENTERY = WIDTH/2

BLACK = (  0,   0,   0)

size = [WIDTH, HEIGHT]
screen = pygame.display.set_mode(size)

# Network
stop = 100;
epsilon = 2000;
fun = 0;
num_of_agents = 4
list_of_agents = []
list_of_edges = [] # For the coordination on the circle

for i in range(num_of_agents-1):
  list_of_edges.append((i,i+1))    # Daisy chain network topology

desired_inter_vehicle_theta = (2*np.pi/num_of_agents)*np.ones((len(list_of_edges),1))
theta_vehicles = np.zeros((num_of_agents,1))
k_coord = 40

# Incidence matrix
B = np.zeros((num_of_agents, len(list_of_edges)))
for idx,edge in enumerate(list_of_edges):
    B[edge[0],idx] =  1
    B[edge[1],idx] = -1

B_dir = np.zeros((num_of_agents, len(list_of_edges)))
for idx,edge in enumerate(list_of_edges):
    B_dir[edge[0],idx] =  1
    B_dir[edge[1],idx] =  0

# Directed vs Undirected graph. They have different convergence properties
# If you want to control with an undirected graph, then uncomment the following.
B_dir = B

for i in range(num_of_agents):
    theta_o = np.pi - (2*np.pi)*np.random.rand(1);
    list_of_agents.append(ag.AgentUnicycle(list(pygame.colordict.THECOLORS.items())[np.random.randint(0,657)][1], i, 1000*np.random.rand(2,1), 40*np.array([[np.cos(theta_o[0])],[np.sin(theta_o[0])]])))

for agent in list_of_agents:
    agent.traj_draw = True

# GVF
ke_circle = 5e-5
kd_circle = 60

xo = 500 # Circle's center (a su vez tengo ya el centro definido lo paso desde fuera).
yo = 500
ro = 50 # radius (esta es mi D la puedo pasar desde fuera).

ck = np.array([xo,yo])

# run simulation
pygame.init()
clock = pygame.time.Clock()
fps = 50
dt = 1.0/fps
time = 0

runsim = True
while(runsim):
    screen.fill(BLACK)

    us = 0 # We keep constant velocity

    # Coordination algorithm on the circle (calculated in a compact way)
    for idx,agent in enumerate(list_of_agents):
        theta_vehicles[idx] = np.arctan2(agent.pos[1]-yo, agent.pos[0]-xo)

    inter_theta = B.transpose().dot(theta_vehicles)
    error_theta = inter_theta - desired_inter_vehicle_theta

    if np.size(error_theta) > 1:
      for i in range(0, np.size(error_theta)):
        if error_theta[i] > np.pi:
          error_theta[i] = error_theta[i] - 2*np.pi
        elif error_theta[i] <= -np.pi:
          error_theta[i] = error_theta[i] + 2*np.pi
    else:
        if error_theta > np.pi:
          error_theta = error_theta - 2*np.pi
        elif error_theta <= -np.pi:
          error_theta = error_theta + 2*np.pi

    dr = -k_coord*B_dir.dot(error_theta)
    # print(la.norm(error_theta))

    # Algorithm gradient estimation. 
    positions_agents = np.zeros((2*num_of_agents,1))
    for idx,agent in enumerate(list_of_agents):
        positions_agents[2*idx] =  np.array([agent.pos[0]]) # Tengo que pasar la posicion, xo, yo y D, y con todo eso los organizo desde el otro lado.
        positions_agents[2*idx + 1] =  np.array([agent.pos[1]]) # Tengo que pasar la posicion, xo, yo y D, y con todo eso los organizo desde el otro lado.
    
    # print(positions_agents)
    # # Seria conveniente hacer que avance aca, porque creo que en caso contrario no me hace la animacion
    # # Voy a tener que separar los codigos, uno para calcular el gradiente 
    if (la.norm(error_theta) < 0.3 and fun < 0.9999):
        gradestfin=grad.computegradient(xo,yo,positions_agents,ro) # Gradiente.
        fun = grad.function(xo,yo) # Gradiente.
        # positions_agents = positions_agents + epsilon*gradestfin
        ck = ck + epsilon*gradestfin # Avance.
        xo = ck[0]
        yo = ck[1]
        stop = sum(abs(gradestfin))
        print("soy fun: ",fun)
        print("gradestfin: ",gradestfin)
        print("soy stop: ",stop)
        print("soy ck: ",ck)
            
    # Guiding vector field
    for idx,agent in enumerate(list_of_agents):
        agent.draw(screen)
        circle_path = gvf.Path_gvf_circle(xo, yo, ro+dr[idx])
        ut = gvf.gvf_control_2D_unicycle(agent.pos, agent.vel, ke_circle, kd_circle, circle_path, 1)
        agent.step_dt(us, ut, dt)

    clock.tick(fps)
    pygame.display.flip()

    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            endtime = pygame.time.get_ticks()
            pygame.quit()
            runsim = False

# Postprocessing
fig = pl.figure(0)
ax = fig.add_subplot(111)
for agent in list_of_agents:
    lp.plot_position(ax, agent)
    ax.axis("equal")
    ax.grid
#circle = pl.Circle((xo, yo), ro, color='r')
#ax.add_patch(circle)

pl.show()
