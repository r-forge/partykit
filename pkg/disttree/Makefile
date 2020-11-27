

vdisttree:=$(shell cat DESCRIPTION | egrep Version | sed 's/Version://g')


.PHONY: codecov check build

codecov: SHELL:=/bin/bash
codecov: 
	Rscript -e "library('covr'); x <- package_coverage(); saveRDS(x, file = '_coverage.rds')"

check: SHELL:=/bin/bash
check:
	-mkdir ../tmp_disttree_inst/ 
	-find ./inst/ -maxdepth 1 -mindepth 1 -not -name tinytest \
    -exec mv '{}' ../tmp_disttree_inst/ \; 
	-cd .. && \
	R CMD build --no-build-vignettes disttree && \
	R CMD check disttree_$(shell printf "%s"${vdisttree}).tar.gz
	-cd disttree
	-mv -vi ../tmp_disttree_inst/* ./inst/
	-rmdir ../tmp_disttree_inst


build: SHELL:=/bin/bash
build:
	-mkdir ../tmp_disttree_inst/ 
	-find ./inst/ -maxdepth 1 -mindepth 1 -not -name tinytest \
    -exec mv '{}' ../tmp_disttree_inst/ \; 
	-cd .. && \
	R CMD build --no-build-vignettes disttree && \
	R CMD INSTALL disttree_$(shell printf "%s"${vdisttree}).tar.gz
	-cd disttree
	-mv -vi ../tmp_disttree_inst/* ./inst/
	-rmdir ../tmp_disttree_inst
