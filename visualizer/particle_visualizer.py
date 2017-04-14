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
	# print("HELLO!!!!")
	args = sys.argv
	if len(args) != 4:
		print("Usage: particle_visualizer.py path/to/particles/csv world_x_bound world_y_bound")
		sys.exit()
	else:
		csv_path = args[1]
		x_bound = int(args[2])
		y_bound = int(args[3])

		robot_radius = float(x_bound)/100.0
		particle_radius = int(x_bound)/200.0

		with open(csv_path, 'r') as csvfile:
			reader = csv.reader(csvfile, delimiter=',')
			# row_counting_reader = csv.reader(csvfile, delimiter=',')
			# row_count = sum(1 for row in row_counting_reader)
			count = 0
			for idx, row in enumerate(reader):
				# print(count)
				# count += 1
				# if count < 900:
				# 	continue
				print(idx)
				pv = ParticleViewer(x_bound, y_bound)
				if len(row) < 2:
					print("No robot provided in line {0} of csv".format(idx))
					sys.exit()
				robot_x = row[0]
				robot_y = row[1]
				pv.add_robot(robot_x, robot_y, radius=robot_radius)
				el_num = 2
				print("Length of row is {0}".format(len(row)))
				while (el_num < len(row) - 1):
					# print(el_num)
					particle_x = row[el_num]
					particle_y = row[el_num + 1]
					pv.add_particle(particle_x, particle_y, radius=particle_radius)
					el_num += 2
				pv.show()
				print("Press enter to continue simulation")
				raw_input()

