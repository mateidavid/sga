SHELL = /bin/bash

MAC_DIR = ${HOME}/code/sga/build/gcc-test/Dev
TAGS = cat-step
OPTS = -s

INPUT_DIR = ..
INPUT_FILES = $(wildcard ${INPUT_DIR}/*.asqg.gz ${INPUT_DIR}/*.asqg)

output_file = $(patsubst %.asqg,%,$(patsubst %.gz,%,$(notdir $(1)))).${TAGS}.mac
TARGETS = $(foreach f,${INPUT_FILES},$(call output_file,${f}))

.PHONY: all list clean

all: ${TARGETS}

list:
	@echo ${TARGETS}

clean:
	@rm -rf ${TARGETS}
#	@rm -rf $(foreach f,${TARGETS},${f} ${f}.log ${f}.stats)

define add_run
$(2): $(1)
	tyme ${MAC_DIR}/mac -d info -c 1000 ${OPTS} -i $(1) -S $$@ --stats-file $$@.stats 2>$$@.log
endef

$(foreach f,${INPUT_FILES},$(eval $(call add_run,${f},$(call output_file,${f}))))
