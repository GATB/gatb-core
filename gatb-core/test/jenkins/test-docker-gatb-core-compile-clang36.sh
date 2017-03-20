# == GATB Compiler machine
# must exist as a docker container on the VM!
COMP_MACHINE=gatb_compiler_clang36

# == enter working directory
DK_WORK_DIR=/builds/workspace/${JOB_NAME}
cd ${DK_WORK_DIR}

# == we have a dedicated directory per BRANCH_TO_BUILD
[ ! -d ${BRANCH_TO_BUILD} ] && { mkdir ${BRANCH_TO_BUILD}; }
cd ${BRANCH_TO_BUILD}

# == we do the git clone (docker container cannot do that, for now)
[ ! -d gatb-core ] && { git clone git+ssh://gatb-ci@scm.gforge.inria.fr/gitroot/gatb-core/gatb-core.git; }

# == we get the appropriate branch to build
cd gatb-core
git checkout ${BRANCH_TO_BUILD}
git pull
cd ..

# == important notice
#    do not delete build directory: it is done by 
#    gatb-compile.sh script, below. In addition
#    user 'ci' won't have permission to do that:
#    build dir being created from the container
#    perspective, it is own by root.

# == we set some variables to prepare volume mount points between container and host
# on host, gatb is here:
DK_MOUNT=${DK_WORK_DIR}/${BRANCH_TO_BUILD}
# from the container, we access source code here:
G_CODE=/tmp/gatb-core-code
# from the container, we prepare build here:
G_BUILD=/tmp/gatb-core-build

# == we run docker to *COMPILE* GATB-Core
docker run --rm --name ${COMP_MACHINE} -e "GIT_PROVIDER=ci" -v ${DK_MOUNT}:${G_CODE} -v ${DK_MOUNT}:${G_BUILD} ${COMP_MACHINE} gatb-compile.sh

# == we run docker to *TEST* GATB-Core
docker run --rm --name ${COMP_MACHINE} -e "GIT_PROVIDER=ci" -v ${DK_MOUNT}:${G_CODE} -v ${DK_MOUNT}:${G_BUILD} ${COMP_MACHINE} gatb-test.sh

