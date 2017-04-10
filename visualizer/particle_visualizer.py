from viewer import ParticleViewer, WeightsViewer
import sys
import pdb
import csv

"""
Input CSV formatting:
Each row is a step:
 - First two value in a row represent the x and y coordinates of the robot.
 - Each two values thereafter represent the x and y coordinates of the particle.

"""

if __name__ == '__main__':
	args = sys.argv
	if len(args) != 4:
		print("Usage: particle_visualizer.py path/to/particles/csv world_x_bound world_y_bound")
		sys.exit()
	else:
		csv_path = args[1]
		x_bound = args[2]
		y_bound = args[3]

	robot_radius = float(x_bound)/100.0
	particle_radius = int(x_bound)/200.0

	with open(csv_path, 'rb') as csvfile:
		reader = csv.reader(csvfile, delimiter=',')
		for idx, row in enumerate(reader):
			pv = ParticleViewer()
			if len(row) < 2:
				print("No robot provided in line {0} of csv".format(idx))
				sys.exit()
			robot_x = row[0]
			robot_y = row[1]
			pv.add_robot(robot_x, robot_y, radius=robot_radius)
			el_num = 2
			while (el_num < len(row)):
				particle_x = row[el_num]
				particle_y = row[el_num + 1]
				pv.add_particle(particle_x, particle_y, radius=particle_radius)
				el_num += 2
			pv.show()
			print("Press enter to continue simulation")
			raw_input()

