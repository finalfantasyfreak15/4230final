CC=/usr/local/apps/cuda/3.2/cuda/bin/nvcc
INCLUDE=-I/usr/local/apps/cuda/3.2/cuda/include \
        -I/usr/local/apps/cuda/SDK2/C/common/inc

#export LD_LIBRARY_PATH=/usr/local/apps/cuda/3.2/cuda/lib64:$LD_LIBRARY_PATH
LIBDIR=-L/usr/local/apps/cuda/SDK2/C/lib 
LIBS=-lcutil

SOURCE=clean_fMRI.cu
EXECUTABLE=clean_fMRI

$(EXECUTABLE): $(SOURCE)
	$(CC) $(INCLUDE) $(LIBDIR) $< -o $@ $(LIBS)

clean:
