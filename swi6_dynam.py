import os
import numpy as np 
from vpython import *
import math
import random

# Diffusivity Constant
D_trans = 0.09#0.0000005
D_rot = 3
kT = 2.5

# Periodic Boundary Condition: 
L = 500
class Particle(sphere): 
	
	def giveNum(self,num):
		self.num = num

	def giveMass(self,mass):
		self.mass = mass

	def getContactPoints(self,theta,phi, theta_offset):

		### Spherical coordinates of points using radius of particle; then convert to Cartesian Coordinates
		self.x1 = self.pos.x + self.radius * cos(theta) * sin(phi)
		self.y1 = self.pos.y + self.radius * sin(theta) * sin(phi)
		self.z1 = self.pos.z + self.radius * cos(phi)

		self.x2 = self.pos.x + self.radius * cos(theta + theta_offset) * sin(phi)
		self.y2 = self.pos.y + self.radius * sin(theta + theta_offset) * sin(phi)
		self.z2 = self.pos.z + self.radius * cos(phi)

class Fuse(compound):


	def giveProperties(self, particle, num, mass):
		self.radius = particle.radius
		self.hasOligomerized = False
		self.numPartners = 0
		self.num = num
		self.mass = mass


	def saveContactPoints(self, cp1,cp2):
		self.cp1 = cp1
		self.cp2 = cp2
	
	def update(self, particle_array = [], delta_t = 1): 


		"""Translational Diffusion """
		#Variance is set to the mean squared displacement 
		delta_x, delta_y, delta_z = np.random.normal(0, 2*D_trans*delta_t, 3)

		force_sum_x, force_sum_y, force_sum_z, force_total = self.calculateEnergy(particle_array)

		self.pos.x = self.pos.x + delta_x + (D_trans/kT *force_sum_x*delta_t)
		self.pos.y = self.pos.y + delta_y + (D_trans/kT *force_sum_y*delta_t)
		self.pos.z = self.pos.z + delta_z + (D_trans/kT *force_sum_z*delta_t)

		""" Rotational Diffusion"""
		# Variance is set to the rotational mean squared displacement
		angle = np.random.normal(0,2*D_rot*delta_t,1) + (D_rot/kT)*self.radius*force_total*delta_t
		self.rotate(angle)

	def calculateEnergy(self, particle_array):
		sigma = 2 * self.radius #optimal distance
		epsilon = 64 # well depth

		hasOligomerized = False
		distance_array_x = []
		distance_array_y = []
		distance_array_z = []
		distance_array_r = []
		particle_num_array = []
		
		""" Get Distances """
		for particle_index in range(0,len(particle_array)):
			if particle_index != particle_array[particle_index].num:
				distance_x, distance_y, distance_z, distance_r = self.calculateDistance(particle_array[particle_index])
				distance_array_x.append(distance_x)
				distance_array_y.append(distance_y)
				distance_array_z.append(distance_z)
				distance_array_r.append(distance_r)
				particle_num_array.append(particle_array[particle_index].num)

		""" Potential Energy of Translation"""
		force_sum_x = 0
		force_sum_y = 0 
		force_sum_z = 0
		force = 0
		for distance_idx in range(0,len(distance_array_x)): 
			if distance_array_r[distance_idx] != 0:
				force = 4*epsilon*((sigma**12/distance_array_r[distance_idx]**13) - (sigma**6/distance_array_r[distance_idx]**7))
				#print(force)
				force_x = (distance_array_x[distance_idx]/distance_array_r[distance_idx]) * force
				force_y = (distance_array_y[distance_idx]/distance_array_r[distance_idx]) * force
				force_z = (distance_array_z[distance_idx]/distance_array_r[distance_idx]) * force

				force_sum_x += force_x
				force_sum_y += force_y
				force_sum_z += force_z

		return force_sum_x, force_sum_y, force_sum_z, force


	def calculateDistance(self, particle2):
		return abs(self.pos.x - particle2.pos.x), abs(self.pos.y - particle2.pos.y), abs(self.pos.z - particle2.pos.z), math.sqrt(((self.pos.x - particle2.pos.x)**2) +((self.pos.y - particle2.pos.y)**2) + (self.pos.z - particle2.pos.z)**2)

class Oligomer(Fuse):

	def giveProperties(self,particle_array_length): 
		self.num = particle_array_length - 1
		self.hasOligomerized = False
		self.radius = 3
		self.numPartners = 0
	def giveMass(self, particle1, particle2):
		self.mass = particle1.mass + particle2.mass
	
	def giveRadius(self, particle1, particle2):
		self.radius = particle1.radius + particle2.radius


	def addPartner(self, particle1, particle2):
		self.numPartners = particle1.numPartners + particle2.numPartners

""" Check Distance Between Contact Points and Spits out Probability that the particles will oligomerize"""
def contactProb(particle1, particle2):
	min_dist = min(math.sqrt((particle1.cp1.pos.x - particle2.cp1.pos.x)**2 + (particle1.cp1.pos.y - particle2.cp1.pos.y)**2 + (particle1.cp1.pos.z - particle2.cp1.pos.z)**2), 
		math.sqrt((particle1.cp1.pos.x - particle2.cp2.pos.x)**2 + (particle1.cp1.pos.y - particle2.cp2.pos.y)**2 + (particle1.cp1.pos.z - particle2.cp2.pos.z)**2),
				math.sqrt((particle1.cp2.pos.x - particle2.cp1.pos.x)**2 + (particle1.cp2.pos.y - particle2.cp1.pos.y)**2 + (particle1.cp2.pos.z - particle2.cp1.pos.z)**2), 
				math.sqrt((particle1.cp2.pos.x - particle2.cp2.pos.x)**2 + (particle1.cp2.pos.y - particle2.cp2.pos.y)**2 + (particle1.cp2.pos.z - particle2.cp1.pos.z)**2))
	return abs(min_dist - 10^-1)/(10^-1)

"""" Check for Contacts before Doing Anything Else""" 
def oligoCheck(particle_array):
	""" Find out if dimerization should occur """
	sigma = 2* particle_array[0].radius #optimal distance
	epsilon = -50 # well depth

	""" Get Distances """
	oligomer_pair_list = []
	for particle1_idx in range(0,len(particle_array)):
		distance_array_x = []
		distance_array_y = []
		distance_array_z = []
		distance_array_r = []
		particle_num_array = []
		for particle2_idx in range(0,len(particle_array)):
			if particle1_idx != particle2_idx:
				distance_x, distance_y, distance_z, distance_r = particle_array[particle1_idx].calculateDistance(particle_array[particle2_idx])
				distance_array_x.append(distance_x)
				distance_array_y.append(distance_y)
				distance_array_z.append(distance_z)
				distance_array_r.append(distance_r)
				particle_num_array.append(particle_array[particle2_idx].num)

		for j in range(0, len(distance_array_r)):
			if distance_array_r[j] <= sigma: 
				oligomer_pair_list.append((particle_array[particle1_idx].num,particle_num_array[j]))


	keep_list = []

	### Account for double counting
	for pair in  oligomer_pair_list:
		if (pair[1],pair[0]) not in keep_list:
			keep_list.append(pair)

	for pair in keep_list: 
		#if particle_array[pair[0]].numPartners > -1:
		#if not particle_array[pair[1]].hasOligomerized:
		temp1_oligomer = Oligomer([particle_array[pair[0]],particle_array[pair[1]]])
		temp2_oligomer = Oligomer([particle_array[pair[0]],particle_array[pair[1]]])
		
		temp1_oligomer.giveProperties(len(particle_array))
		temp2_oligomer.giveProperties(len(particle_array))
		temp1_oligomer.addPartner(particle_array[pair[0]],particle_array[pair[1]])
		temp1_oligomer.addPartner(particle_array[pair[0]],particle_array[pair[1]])

		particle_array[pair[0]] = temp1_oligomer
		particle_array[pair[1]] = temp2_oligomer
		particle_array[pair[0]].hasOligomerized = True

	keep_list = [particle for particle in particle_array if particle.hasOligomerized == False]

	for particle in particle_array:
		if particle.hasOligomerized:
			particle.visible = False
			del particle

	return keep_list


"""Main Simulation - The Conductor! """
particle_array = []
for i in range(0,20):
	#temp_particle = Particle(pos = vector(i,0,0), radius = 0.8, color = color.white)

	temp_particle = Particle(pos = vector(random.uniform(-10,10),random.uniform(-10,10),random.uniform(-10,10)), radius = 0.8, color = color.white)
	temp_particle.getContactPoints(theta = 0, phi = 30, theta_offset = 20)
	contactPoint1 = Particle(pos = vector(temp_particle.x1,temp_particle.y1,temp_particle.z1), radius = 0.2, color = color.green)
	contactPoint2 = Particle(pos = vector(temp_particle.x2,temp_particle.y2,temp_particle.z2), radius = 0.2, color = color.blue)
	
	temp_particle_fusion = Fuse([temp_particle, contactPoint1])
	temp_particle_fusion = Fuse([temp_particle_fusion, contactPoint2])
	temp_particle_fusion.saveContactPoints(contactPoint1, contactPoint2)
	temp_particle_fusion.giveProperties(temp_particle, num =i, mass =2)

	particle_array.append(temp_particle_fusion)

for __ in range(0,1000):
	rate(400) 
	print("Now Checking for Oligomers")
	particle_array = oligoCheck(particle_array)
	print(particle_array)
	print("Oligomers made!")
	for particle in particle_array:
		particle.update(particle_array)
	print("New Move!")
print("Done")
#input("next ...")