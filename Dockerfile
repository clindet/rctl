FROM bioinstaller/bioinstaller

## This handle reaches Jianfeng
MAINTAINER "Jianfeng Li" lee_jianfeng@life2cloud.com

Add  . /tmp/ngsjs

RUN su opencpu \
    && /bin/bash \
    && export PATH=/home/opencpu/opt/miniconda2/bin:$PATH \
    && conda install go nodejs \
    && conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free \
    && conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge \
    && conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/msys2 \
    && echo 'export NODE_PATH="/home/opencpu/opt/miniconda2/lib/node_modules/"\n' >> /etc/profile \
    && npm install -g npm \
    && npm install -g yarn \
    && echo 'export PATH=$PATH:~/.yarn/bin/\n' >> /etc/profile \
    && yarn global add /tmp/ngsjs \
    && rm -rf /tmp/ngsjs\
    && rdeps

CMD runuser -s /bin/bash -l opencpu -c '. /etc/profile; ngsjs'
