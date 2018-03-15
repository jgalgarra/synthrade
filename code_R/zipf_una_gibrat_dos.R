exporters <- seq(1,100)
importers <- seq(1,100)
alfa <- 1.4
beta <- 1.5
P_exporters <- 1/(exporters^alfa)
P_exporters <- P_exporters/sum(P_exporters)

P_importers <- 1/(importers^alfa)
P_importers <- P_importers/(importers^beta)

P_tot <- P_importers %o% P_exporters