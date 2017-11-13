#include "encoder.h"
#include <Eigen/Eigen>	//Includes Dense and Sparse header files (the whole Eigen library)
#include <math.h>
#include <iostream>
#include <algorithm>

//-------------------------------------------- Function definitions ---------------------------------------------

/*Construct an octree of depth b, based on sorted (ascending order) b-bit Morton codes of the input point cloud.
This octree will have b + 1 levels in total, from level 0 (root) to level b (leaves).*/
void constructOctree(const int b, octree& my_OT)
{


	return;
}

/*For each occupied octree cell at each octree level, compute the x, y, z coordinates of each corner of this cell,
then extract only the unique corners at each level and the pointer to the associated Bezier control point for
each unique corner.*/
void computeCornerCoordinates(const int b, const int max_lvl, corners& input_corners, octree& my_OT)
{
	for (int lvl = 0; lvl < max_lvl; lvl++)
	{
		//Find the (x, y, z) coordinates of the origin of each occupied octree block at the current level (consider the origin corner 1)
		Eigen::Matrix<double, my_OT.NodeCount(lvl), 3> corners1 = ;

		//Replicate each row of corners1 7 times (so that the resulting matrix can be added directly to offsets_from_origin_rep, below)
		Eigen::Matrix<double, 7*corners1.rows(), 3> corners1_rep = ;

		//Find the (x, y, z) coordinates of the other 7 corners of each occupied octree block at the current level
		Eigen::Matrix<double, 8, 3> offsets_from_origin;
		offsets_from_origin << pow(2, (b + 1 - lvl)), 0, 0,
			pow(2, (b + 1 - lvl)), pow(2, (b + 1 - lvl)), 0,
			0, pow(2, (b + 1 - lvl)), 0,
			0, 0, pow(2, (b + 1 - lvl)),
			pow(2, (b + 1 - lvl)), 0, pow(2, (b + 1 - lvl)),
			pow(2, (b + 1 - lvl)), pow(2, (b + 1 - lvl)), pow(2, (b + 1 - lvl)),
			0, pow(2, (b + 1 - lvl)), pow(2, (b + 1 - lvl));
		Eigen::Matrix<double, corners1_rep.rows(), 3> offsets_from_origin_rep = ;
		Eigen::Matrix<double, offsets_from_origin_rep.rows(), 3> corners2_8 = corners1_rep + offsets_from_origin_rep;	//Corners 2-8 for each octree cell

		/*Reorganize all the corner coordinates for the current octree level, so that all 8 corners for one octree
		cell have their coordinates listed one after another (i.e, in an 8 x 3 matrix), inside corner_coords.*/
		Eigen::Matrix<double, (corners1.rows() + corners2_8.rows()), 3> corner_coords = ;

		//Extract the coordinates of only the unique corners at the current octree level
		for (int r = 0; r < corner_coords.rows(); r++)
		{
			//Get the r-th row of coordinates
			Eigen::Array<double, 1, 3> current_row = corner_coords.data.row(r);

			//If the current row has not already been stored inside unique_coords, store it (otherwise ignore it)
			int* p = std::find(input_corners.unique_coords(lvl).data.row(0), input_corners.unique_coords(lvl).data.row(), current_row);
			if (p != input_corners.unique_coords(lvl).row())
			{
				input_corners.unique_coords(lvl)() = current_row;
			}
		}
		

		//For each row of coordinates in corner_coords, store a pointer to the corresponding row in unique_coords

	}

	return;
}

/*For each occupied octree cell at each octree level, extract the occupied voxel coordinates associated with that
cell, as well as the corresponding normals and centroids if the input point cloud has those.*/
extractOccupiedVoxels()
{

}

/*For each unique corner of the occupied blocks at each octree level, compute the Bezier control point (SDF
value) for this corner.*/
void computeControlPoints(const int ctrlpt_thresh, const int max_lvl, const int b, corners& input_corners, Eigen::Matrix<Eigen::VectorXd, int depth, 1>& control_points)
{
	/*The number of control points inside each vector in control_points will be equal to the number of 
	unique corners at the corresponding octree level*/
	for (int i = 0; i < (b + 1); i++)
	{
		control_points(i).resize(unique_coords(i).rows(), 1);
	}

	for (int lvl = 0; lvl < max_lvl; lvl++)
	{
		if (lvl < thresh)
		{
			/*For each unique corner of the occupied blocks at the current octree level, find the nearest 
			occupied voxel in any of the occupied octree cells at the current level that share this corner*/
			for (int c = 0; c < unique_coords(lvl).rows(); c++)
			{

			}	//End c

		}
		else if (lvl >= thresh)
		{
			//If we are at the voxel level
			if (lvl == b + 1)
			{

			}
			//If we are not at the voxel level
			else
			{

			}
		}
	}

	return;
}

/*For each octree level after start_lvl, compute the wavelet coefficients associated with the unique corners of
the occupied blocks at these levels, and compute the reconstructed control points at these levels that result
after wavelet coefficient quantization and dequantization (done to prevent quantization error accumulating).*/
waveletAnalysis()
{

}

/*Going bottom-up from the voxel level, prune branches (occupancy codes) of octree cells that contain all zero
(or within the chosen threshold) wavelet coefficients throughout the branch. This will result in leaf cells at
variable octree levels.*/
pruneOctree()
{

}

/*Prune the wavelet coefficient tree to match the pruned octree: wavelet coefficients will then exist only for
the unpruned occupied octree cells at each level.*/
pruneWaveletCoeffTree()
{

}

/*Encode the data that will be transmitted to the decoder, and write it to a bitstream.*/
encode()
{

}

//------------------------------------------ Class method definitions -------------------------------------------



//------------------------------------------- Begin program execution -------------------------------------------

int main(int argc, char* argv[])
{
	std::string ptcloud_file;	//Full path to input point cloud
	int b;						//Octree depth (root = level 0; leaves = level b)
	int start_lvl;				//Octree level for which to transmit control points to the decoder, and from which to start wavelet analysis
	int max_lvl;				//Highest octree level for which Bezier control points will be computed at the encoder
	int q_stepsize;				//Quantization stepsize for uniform scalar quantization of control points at start_lvl and all wavelet coefficients
	int ctrlpt_thresh;			//Threshold above which control points will be computed as average dot products for each unique corner
	int zero_threshold;			//Threshold at/below which wavelet coefficients are considered zero (this information is required for pruning)

	//Check for correct number of input arguments to the command line (need 8 input arguments, the first one being the executable name)
	if (argc < 8)
	{
		std::cerr << "Usage: " << argv[0] << " ptcloud_file b start_lvl max_lvl q_stepsize ctrlpt_thresh zero_threshold" << std::endl;

		return 1;
	}
	else
	{
		ptcloud_file = argv[1];
		b = std::stoi(argv[2]);
		start_lvl = std::stoi(argv[3]);
		max_lvl = std::stoi(argv[4]);
		q_stepsize = std::stoi(argv[5]);
		ctrlpt_thresh = std::stoi(argv[6]);
		zero_threshold = std::stoi(argv[7]);
	}

	//Read in input point cloud (SVO format)
	

	/*Construct an octree of depth b, based on sorted (ascending order) b-bit Morton codes of the input point cloud.
	This octree will have b + 1 levels in total, from level 0 (root) to level b (leaves).*/
	octree my_OT;
	constructOctree(b, my_OT);

	/*For each occupied octree cell at each octree level, compute the x, y, z coordinates of each corner of this cell,
	then extract only the unique corners at each level and the pointer to the associated Bezier control point for
	each unique corner.*/
	corners input_corners;
	computeCornerCoordinates(b, max_lvl, input_corners, my_OT);

	/*For each occupied octree cell at each octree level, extract the occupied voxel coordinates associated with that
	cell, as well as the corresponding normals and centroids if the input point cloud has those.*/
	extractOccupiedVoxels();

	/*For each unique corner of the occupied blocks at each octree level, compute the Bezier control point (SDF
	value) for this corner.*/
	Eigen::Matrix<Eigen::VectorXd, (b + 1), 1> control_points;		//Initialize a matrix of (b + 1) vectors (of doubles), to store the control points at each octree level
	computeControlPoints(ctrlpt_thresh, max_lvl, &input_corners, &control_points);

	/*For each octree level after start_lvl, compute the wavelet coefficients associated with the unique corners of
	the occupied blocks at these levels, and compute the reconstructed control points at these levels that result
	after wavelet coefficient quantization and dequantization (done to prevent quantization error accumulating).*/
	waveletAnalysis();

	/*Going bottom-up from the voxel level, prune branches (occupancy codes) of octree cells that contain all zero
	(or within the chosen threshold) wavelet coefficients throughout the branch. This will result in leaf cells at
	variable octree levels.*/
	pruneOctree();

	/*Prune the wavelet coefficient tree to match the pruned octree: wavelet coefficients will then exist only for
	the unpruned occupied octree cells at each level.*/
	pruneWaveletCoeffTree();

	/*Encode the data that will be transmitted to the decoder, and write it to a bitstream.*/
	encode();

	return 0;
}