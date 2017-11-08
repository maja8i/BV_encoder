//------------------------------------ Function definitions -------------------------------------

/*Construct an octree of depth b, based on sorted (ascending order) b-bit Morton codes.
This octree will have b + 1 levels in total, from level 0 (root) to level b (leaves).*/
construct_octree()
{

}

/*For each occupied octree cell at each octree level, compute the x, y, z coordinates of each
corner of this cell, then extract only the unique corners at each level.*/
compute_corner_coordinates()
{

}

/*For each occupied octree cell at each octree level, extract the occupied voxel coordinates
associated with that cell, as well as the corresponding normals and centroids if the input
point cloud has those.*/
extract_occupied_voxels()
{

}

/*For each unique corner of the occupied blocks at each octree level, compute the Bezier
control point (SDF value) for this corner.*/
compute_control_points()
{

}

/*For each octree level after start_lvl, compute the wavelet coefficients associated with
the unique corners of the occupied blocks at these levels, and compute the reconstructed
control points at these levels that result after wavelet coefficient quantization and
dequantization (done to prevent quantization error from accumulating).*/
wavelet_analysis()
{

}

/*Going bottom-up from the voxel level, prune branches (occupancy codes) of octree cells
that contain all zero (or within the chosen threshold) wavelet coefficients throughout
the branch. This will result in leaf cells at variable octree levels.*/
prune_octree()
{

}

/*Prune the wavelet coefficient tree to match the pruned octree: wavelet coefficients will
then exist only for the unpruned occupied octree cells at each level.*/
prune_wavelet_coeff_tree()
{

}

/*Encode the data that will be transmitted to the decoder, and write it to a bitstream.*/
encode()
{

}

//---------------------------------- Class method definitions -----------------------------------



//----------------------------------- Begin program execution -----------------------------------

int main(void)
{
	construct_octree();

	compute_corner_coordinates();

	extract_occupied_voxels();

	compute_control_points();

	wavelet_analysis();

	prune_octree();

	prune_wavelet_coeff_tree();

	encode();

	return 0;
}