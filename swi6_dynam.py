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

	def giveRadius(self,particle):
		self.radius = particle.radius

	def giveNum(self,num):
		self.num = num
	
	def giveMass(self,mass):
		self.mass = mass
	
	def calculateEnergy(self, particle_array):
		sigma = 64* self.radius #optimal distance
		epsilon = -50 # well depth

		distance_array = []


		for particle_index in range(0,length(particle_array)):
			if particle_index != self.num:
				distances_array.append(calculateDistance(self, particle))
		
		""" Potential Energy of Translation"""
		force_sum = 0
		for distance in distances_array: 
			force_sum += 4*epsilon*((sigma/distance)**12 - (sigma/distance)**6)

		return force_sum 


	def update(self, particle_array = [], delta_t = 1): 


		"""Translational Diffusion """
		#Variance is set to the mean squared displacement 
		delta_x, delta_y, delta_z = np.random.normal(0, 2*D_trans*delta_t, 3)

		force_sum_x, force_sum_y, force_sum_z, force_total = self.calculateEnergy(particle_array)

		self.pos.x = self.pos.x + delta_x + (D_trans/kT *force_sum_x*delta_t)
		self.pos.y = self.pos.y + delta_y + (D_trans/kT *force_sum_y*delta_t)
		self.pos.z = self.pos.z + delta_z + (D_trans/kT *force_sum_z*delta_t)



		"""Periodic Boundary Conditions"""
		# if self.pos.x < -1*L: 
		# 	self.pos.x += 0
		# 	#self.pos.y += L
		# 	#self.pos.z += L
		# if self.pos.y < -1*L:  
		# 	#self.pos.x += L
		# 	self.pos.y += 0
		# 	#self.pos.z += L
		# if self.pos.z < -1*L:  
		# 	#self.pos.x += L
		# 	#self.pos.y += L
		# 	self.pos.z += 0

		# if self.pos.x > L: 
		# 	self.pos.x -= 0
		# 	#self.pos.y -= L
		# 	#self.pos.z -= L
		# if self.pos.y > L:  
		# 	#self.pos.x -= L
		# 	self.pos.y -= 0
		# 	#self.pos.z -= L
		# if self.pos.z > L:  
		# 	#self.pos.x -= L
		# 	#self.pos.y -= L
		# 	self.pos.z -= 0
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
			if particle_index != self.num:
				distance_x, distance_y, distance_z, distance_r = self.calculateDistance(particle_array[particle_index])
				distance_array_x.append(distance_x)
				distance_array_y.append(distance_y)
				distance_array_z.append(distance_z)
				distance_array_r.append(distance_r)
				particle_num_array.append(particle_array[particle_index].num)
				#print("Distance: ", math.sqrt(distance_x**2 + distance_y**2 + distance_z**2))
				# print(math.sqrt(distance_x**2 + distance_y**2 + distance_z**2) <0.5)
		
		""" Find out if dimerization should occur """
		oligomer_pair_list = []
		for j in range(0, len(distance_array_r)):
			if distance_array_r[j] <= sigma and particle_num_array[j] != self.num: 
				oligomer_pair_list.append((self.num,particle_num_array[j]))
		
		print(particle_array)
		
		for pair in oligomer_pair_list:
			print(pair)
			temp_oligomer = Oligomer([particle_array[pair[0]],particle_array[pair[1]]])
			
			particle_array[pair[0]].visible = False
			particle_array[pair[1]].visible = False

			temp_oligomer.giveMass(particle_array[pair[0]], particle_array[pair[1]])
			temp_oligomer.giveRadius(particle_array[pair[0]], particle_array[pair[1]])
			temp_oligomer.giveNum(particle_array[pair[0]], particle_array[pair[1]])


			particle_array.append(temp_oligomer)
		for pair in oligomer_pair_list:
			del particle_array[pair[0]]
			del particle_array[pair[1] - 1]

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

	def giveNum(self,particle1,particle2): 
		self.num = 2

	def giveMass(self, particle1, particle2):
		self.mass = particle1.mass + particle2.mass
	
	def giveRadius(self, particle1, particle2):
		self.radius = particle1.radius + particle2.radius

"""Main Simulation - The Conductor! """
particle_array = []
for i in range(0,2):
	temp_particle = Particle(pos = vector(random.uniform(-1,2),random.uniform(-1,2),random.uniform(-1,3)), radius = 0.8, color = color.white)
	temp_particle.getContactPoints(theta = 0, phi = 30, theta_offset = 20)
	contactPoint1 = Particle(pos = vector(temp_particle.x1,temp_particle.y1,temp_particle.z1), radius = 0.1, color = color.green)
	contactPoint2 = Particle(pos = vector(temp_particle.x2,temp_particle.y2,temp_particle.z2), radius = 0.1, color = color.blue)
	
	temp_particle_fusion = Fuse([temp_particle,contactPoint1])
	temp_particle_fusion = Fuse([temp_particle_fusion, contactPoint2])
	temp_particle_fusion.giveRadius(temp_particle)
	temp_particle_fusion.giveNum(i)
	temp_particle_fusion.giveMass(2)

	particle_array.append(temp_particle_fusion)
for __ in range(0,2000):
	rate(10) 
	for particle in particle_array:
		print(particle_array)
		particle.update(particle_array)

	print("New Move!")
print("Done")
#input("next ...")