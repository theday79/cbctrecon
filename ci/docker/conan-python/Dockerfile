FROM ubuntu:focal
# focal = 20.04

MAINTAINER Andreas Gravgaard Andersen <andreasga22@gmail.com>

USER root

RUN echo "deb http://us.archive.ubuntu.com/ubuntu/ focal multiverse" >> /etc/apt/sources.list \
    && apt update \
    && echo "Install dependencies for conan" \
    && apt install -y wget sudo git gcc make \
        zlib1g-dev libffi-dev libbz2-dev libssl-dev libsqlite3-dev libreadline-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && echo "Add user with group sudo" \
    && groupadd 1001 -g 1001 \
    && groupadd 1000 -g 1000 \
    && groupadd 2000 -g 2000 \
    && groupadd 999 -g 999 \
    && useradd -ms /bin/bash user -g 1001 -G 1000,2000,999 \
    && printf "user:user" | chpasswd \
    && adduser user sudo \
    && printf "user ALL= NOPASSWD: ALL\\n" >> /etc/sudoers

RUN mkdir /opt/cmake && cd /opt/cmake \
    && wget -nv https://cmake.org/files/v3.17/cmake-3.17.2-Linux-x86_64.sh \
    && sh cmake-3.17.2-Linux-x86_64.sh --prefix=/opt/cmake --skip-license \
    && ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake \
    && ln -s /opt/cmake/bin/ctest /usr/local/bin/ctest \
    && rm cmake-3.17.2-Linux-x86_64.sh \
    && cmake --version

ENV PYENV_ROOT=/opt/pyenv \
    PATH=/opt/pyenv/shims:${PATH}

RUN wget --no-check-certificate --quiet -O /tmp/pyenv-installer https://github.com/pyenv/pyenv-installer/raw/master/bin/pyenv-installer \
    && chmod +x /tmp/pyenv-installer \
    && /tmp/pyenv-installer \
    && rm /tmp/pyenv-installer \
    && update-alternatives --install /usr/bin/pyenv pyenv /opt/pyenv/bin/pyenv 100 \
    && PYTHON_CONFIGURE_OPTS="--enable-shared" pyenv install 3.7.1 \
    && pyenv global 3.7.1 \
    && pip install -q --upgrade --no-cache-dir pip \
    && pip install -q --no-cache-dir conan conan-package-tools \
    && chown -R user:1001 /opt/pyenv \
    # remove all __pycache__ directories created by pyenv
    && find /opt/pyenv -iname __pycache__ -print0 | xargs -0 rm -rf \
    && update-alternatives --install /usr/bin/python python /opt/pyenv/shims/python 100 \
    && update-alternatives --install /usr/bin/python3 python3 /opt/pyenv/shims/python3 100 \
    && update-alternatives --install /usr/bin/pip pip /opt/pyenv/shims/pip 100 \
    && update-alternatives --install /usr/bin/pip3 pip3 /opt/pyenv/shims/pip3 100

USER user

RUN mkdir -p /home/user/.conan \
    && printf 'eval "$(pyenv init -)"\n' >> ~/.bashrc \
    && printf 'eval "$(pyenv virtualenv-init -)"\n' >> ~/.bashrc

