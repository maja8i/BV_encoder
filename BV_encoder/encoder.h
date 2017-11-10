#pragma once

//------------------------------------ Function declarations ------------------------------------

/*Construct an octree of depth b, based on sorted (ascending order) b-bit Morton codes of
the input point cloud. This octree will have b + 1 levels in total, from level 0 (root) to
level b (leaves).*/
void constructOctree(const int b, octree& my_OT);

/*For each occupied octree cell at each octree level, compute the x, y, z coordinates of each
corner of this cell, then extract only the unique corners at each level and the pointer to the
associated Bezier control point for each unique corner.*/
void computeCornerCoordinates(const int b, const int max_lvl, corners& input_corners, octree& my_OT);

/*For each occupied octree cell at each octree level, extract the occupied voxel coordinates
associated with that cell, as well as the corresponding normals and centroids if the input
point cloud has those.*/
extractOccupiedVoxels();

/*For each unique corner of the occupied blocks at each octree level, compute the Bezier
control point (SDF value) for this corner.*/
void computeControlPoints(const int ctrlpt_thresh, const int max_lvl, const int b, corners& input_corners, Eigen::Matrix<Eigen::VectorXd, int depth, 1>& control_points);

/*For each octree level after start_lvl, compute the wavelet coefficients associated with
the unique corners of the occupied blocks at these levels, and compute the reconstructed 
control points at these levels that result after wavelet coefficient quantization and
dequantization (done to prevent quantization error from accumulating).*/
waveletAnalysis();

/*Going bottom-up from the voxel level, prune branches (occupancy codes) of octree cells
that contain all zero (or within the chosen threshold) wavelet coefficients throughout
the branch. This will result in leaf cells at variable octree levels.*/
pruneOctree();

/*Prune the wavelet coefficient tree to match the pruned octree: wavelet coefficients will
then exist only for the unpruned occupied octree cells at each level.*/
pruneWaveletCoeffTree();

/*Encode the data that will be transmitted to the decoder, and write it to a bitstream.*/
encode();

//------------------------------------ Structure definitions ------------------------------------

struct corners
{
	//Matrix of (b + 1) matrices (of doubles), to store only the coordinates of the unique corners at each octree level
	Eigen::Matrix<Eigen::MatrixXd, (b + 1), 1> unique_coords;

	/*Matrix of (b + 1) vectors (of doubles), containing pointers to the Bezier control points that are associated
	with the corners in corner_coords. There will be 1 pointer per corner (8 pointers per occupied octree block) at each 
	octree level.*/
	Eigen::Matrix<Eigen::VectorXd, (b + 1), 1> ctrl_pts_pointers;
};

struct octree
{

};

//-------------------------------------- Class definitions --------------------------------------


