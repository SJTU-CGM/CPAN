export SRC_DIR	= src
SUBDIRS = $(SRC_DIR)/perl $(SRC_DIR)/R

all: 
	@mkdir -p bin;\
	for subdir in $(SUBDIRS); \
	do \
	(cd $$subdir && make); \
	done
.PHONY: all

clean:
	@for subdir in $(SUBDIRS); \
	do \
	(cd $$subdir && make clean); \
	done
