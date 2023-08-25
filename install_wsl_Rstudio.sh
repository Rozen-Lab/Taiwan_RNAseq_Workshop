sudo apt-get update
sudo apt-get upgrade -y
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
sudo apt install -y r-base r-base-core r-recommended r-base-dev gdebi-core build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev libfontconfig1-dev libhdf5-dev
# install Rstudio server version
chmod 644 .wget-hsts #if 777 the wget may not work # error: Will not apply HSTS
sudo apt-get install gdebi-core
wget https://download2.rstudio.org/server/focal/amd64/rstudio-server-2023.06.1-524-amd64.deb
sudo gdebi rstudio-server-2023.06.1-524-amd64.deb
sudo rstudio-server start

# check on the local website
# username and password will be Linux username and password
http://localhost:8787/
