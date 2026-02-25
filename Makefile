# 0000000000000000000000000000000000000000000000000000000000000
# This file is part of TREKIS-4
# available at: https://github.com/N-Medvedev/TREKIS-4
# The makefile calles for another makefile within Source_files directory
# which compiles and makes an executable called TREKIS.x
# This file was written by N.Medvedev 
# in 2020-2026
#----------------------------------------------------- 
# Define a timestamp function
# To pass variables into the next make-file:
export


# Call makefile within the Source_files directory:
subsystem:
	@echo "**********************************************"
	@echo "Starting compilation at"
	@echo $(shell date)
	@echo "**********************************************"

	cd Source_files && $(MAKE)

# Copy created executable into the parent directory:
	scp -r Source_files/*.x .

# Delete executable from the Source_files directory:
	rm -r Source_files/*.x

	@echo "**********************************************"
	@echo "Finished compilation at"
	@echo $(shell date)
	@echo "**********************************************"

# Clean all compiled files:
clean:
	rm -f Source_files/*.o
	rm -f Source_files/*.mod
	rm -f Source_files/*.obj
	rm -f Source_files/*.yaml
	rm -f Source_files/*.optrpt
	rm -f Source_files/*.x

# Clean all results:
cleanresults:
	rm -f OUTPUT_*

# Clean all compiled files and the results:
veryclean:
	rm -f Source_files/*.o
	rm -f Source_files/*.mod
	rm -f Source_files/*.obj
	rm -f Source_files/*.x
	rm -f *.x
	rm -f OUTPUT_*
