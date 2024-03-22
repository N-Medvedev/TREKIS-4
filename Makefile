# The makefile calles for another makefile within Source_files directory
# which compiles and makes an executable called TREKIS.x
# This file was written by N.Medvedev 
# in 2020-2021
#----------------------------------------------------- 
# To pass variables into the next make-file:
export

# Call makefile within the Source_files directory:
subsystem:
	cd Source_files && $(MAKE)

# Copy created executable into the parent directory:
	scp -r Source_files/TREKIS.x TREKIS.x

# Delete executable from the Source_files directory:
	rm -r Source_files/TREKIS.x

clean:
	cd Source_files
	rm -f *.o
	rm -f *.mod
	rm -f TREKIS.x
