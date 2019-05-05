fcm - Zhiyang Huang (2018)

------------------------------------

This code implements the algorithm described in

Huang, Zhi Yang, et al. "Repairing Inconsistent Curve Networks on Non-parallel Cross-sections." Computer Graphics Forum. Vol. 37. No. 2. 2018.

The primary function of this program is to restore the inconsistent cross-sections. Giving input of a set of inconsistent contours, the program builds up the implicit domain (triangulated planes sharing intersection lines) and modifies the implicit representation (sign distance function) for repairing the inconsistency.

Currently, the code is only tested on Mac OS X 10.10 or above, it should work at other platforms with minor modification.


BUILDING
======================================================================================================


The code has only three dependencies: 1)Eigen/Eigen3,   2)Gurobi,  3)Suite-Sparse

1) http://eigen.tuxfamily.org/index.php?title=Main_Page

2) http://www.gurobi.com/

3) http://faculty.cse.tamu.edu/davis/suitesparse.html

After download/install the three dependencies, please modify the relevant path according to your installation in the CMakeList.txt, e.g.

Eigen include file path:    SET(EIGEN_INCLUDE_DIRS "/opt/local/include/")

Gurobi include file and lib path:
SET(GUROBI_INCLUDE_DIRS "/Library/gurobi702/mac64/include/")
SET(GUROBI_LIB_DIRS "/Library/gurobi702/mac64/lib/")

Suite-Sparse lib path:  SET(SUITESPARSE_LIB_DIR "/usr/local/Cellar/suite-sparse/5.2.0/lib/")

Then go to the folder EG2018, build the Cmake file and make:
$cmake .
$make

In the EG2018 directory, there should be an executable called "fcm" (or "fcm.exe" on Windows if it is successfully built).


RUNNING
======================================================================================================

To run the code from the command line, type:

$./fcm -i input_file_name -o output_file_path [-w] [-s line_step] [-l lambda] [-p active_portion] [-m max_iteration] [-t triangle_command]

where:
1. -i: followed by the path of the input file. input_file_name is a path to the input file. currently, support file format includes ".contour". See https://www.cse.wustl.edu/~taoju/lliu/paper/ctr2suf/program.html for the description of the format.

2. -o: followed by the path of the output path. output_file_path is a path to the folder for generating output files.

3. -w: optional argument. If -w is included in the command line, the program output the intermedium results, including the input (mediumresult_input.obj) and the triangulated planes (mediumresult_tri.obj).

4. -s: optional argument. Followed by a float number indicating the sampling steps (e.g., the distance between two sample points) in the intersection lines. Default 0.02, you should set this according to your inputs, you can check the triangulated result in mediumresult_tri.obj by turn on -w option

5. -l: optional argument. Followed by a float number indicating the lambda which balances the energy (see the paper for details). Default 0.001, you should set and tune this number according to your inputs.

6. -p: optional argument. Followed by a float number indicating the portion of active vertices outside those inconsistent vertices, See the Experimental part of the paper for details. Default 0.1.

7. -m: optional argument. Followed by an integer number indicating the maximum number of update in the optimization step. Default 30.

8. -t: optional argument. Followed by a string for calling the Triangle library, which is used for triangulating planes. Default "pzQYYa0.001". See  https://www.cs.cmu.edu/~quake/triangle.html for details. You should modify it according to your input. Please carefully set up the constraints othervise it will fail.

Few examples is placed at data folder for testing:
1. $./fcm -i data/Bermano/Bermano.contour -o data/Bermano/                 -m 30 -w -p 0.1 -s 0.5      -t pzQYYa1         -l 0.005
2. $./fcm -i data/mousebrain/mousebrain.contour -o data/mousebrain/    -m 30 -w -p 0.1 -s 0.02    -t pzQYYa0.001  -l 0.05
3. $./fcm -i data/ferretbrain/ferretbrain.contour -o data/ferretbrain/            -m 30 -w -p 0.1 -s 1        -t pq28zYYa0.8   -l 0.001
4. $./fcm -i data/newliver/newliver.contour -o data/newliver/ 			   -m 30 -w -p 0.1 -s 0.03   -t pqzYYa0.001   -l 0.005





