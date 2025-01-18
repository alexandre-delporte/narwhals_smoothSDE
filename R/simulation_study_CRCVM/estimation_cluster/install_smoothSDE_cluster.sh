source /applis/site/guix-start.sh
export SSL_CERT_DIR=$(guix build nss-certs)/etc/ssl/certs
R -e "library(devtools);install_github('alexandre-delporte/smoothSDE')"

