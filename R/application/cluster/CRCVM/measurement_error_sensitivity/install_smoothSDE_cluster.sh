source /applis/site/guix-start.sh
export SSL_CERT_DIR=$(guix build nss-certs)/etc/ssl/certs
export R_LIBS_USER=~/R/library
mkdir -p ~/R/library
R -e "library(devtools); install_github('alexandre-delporte/smoothSDE', lib='~/R/library')"

