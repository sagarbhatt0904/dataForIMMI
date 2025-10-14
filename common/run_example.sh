NPROCS= <NUM_MPI_RANKS>
ROOTDIR=<PATH_TO>/am_<LOAD_DIR>
SUB=<PATH_TO>/am_<LOAD_DIR>/<INIT_CAVITY_SIZE>
DEER=<PATH_TO>/deer-opt
INPUT=<PATH_TO>/common/y.i
FILE=<PATH_TO>/microstructureData/AMmicro_generic.e
INTPROP=$SUB/gbprop.i
MODELFILE=<PATH_TO>/common/316H_simple.xml
EULER=<PATH_TO>/microstructureData/ori_byBlock.tex


mpirun -n $NPROCS $DEER -i $INTPROP $INPUT Mesh/base/file=$FILE UserObjects/euler_angle_file/prop_file_name=$EULER UserObjects/euler_angle_file/read_type=block UserObjects/euler_angle_file/nblock=3309 Materials/stress/database=$MODELFILE load=130
