# CS302-GJK_algorithm
This project is an implementation of Gilbert-Johnson and Keerthi's algorithm as a part of our CS302-Analysis and Design of algorithm course.
Mainly focused on collision detection of 3D objects.


Compile and run:

	g++ -std=c++17 gjk.cpp
	./a.out

[Inputs](/inputs) folder contains the inputs given to the program. Any input file is expected to be located
in this folder. Input format goes as : 
	line 1 : number of vertices for polyhedra 1 - N1
	N1 lines follow each containing 3 space seperated numbers denoting each vertex.
	line 2+N1: number of vertices for polyhedra 2 - N2
	N2 lines follow each containing 3 space seperated numbers denoting each vertex.

[Blender](/blender) folder contains the .blend files created in blender. We have creted input test cases in
this software. After creating two polyhedra we can import a .obj file which will contain the
vertices of both the polyhedra. We copied this data into .txt files in inputs folder.
Our implementation gjk.cpp takes these vertices as inputs.	
