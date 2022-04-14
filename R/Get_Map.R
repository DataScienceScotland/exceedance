
if (!is.element('sf', installed.packages()[,1]))
  install.packages(p, dep = TRUE)
require('sf', character.only = TRUE)


if(file.exists("Scottish_Map.RData")){
  load("Scottish_Map.RData")
}else{
  M1=read_sf('R/GB_boundaries/Local_Authority_Districts__December_2017__Super_Generalised_Clipped_Boundaries_in_Great_Britain.shp');
  Z=rep(0, length(M1$lad17nm))

  tmp_map_place_names = gsub(" ", "", M1$lad17nm)
  tmp_map_place_names = gsub("-", "", tmp_map_place_names)
  tmp_Place_names = gsub(" ", "", Place_names)
  tmp_Place_names = gsub("-", "", tmp_Place_names)

  MapX = c()
  MapY = c()
  Map2PlaceName = c()
  I = c()
  for (i in 1:length(M1$lad17nm)){
    if (tmp_map_place_names[i] %in%  tmp_Place_names){
      MapX = c(MapX, mean(st_coordinates(M1[i,])[,'X'], na.rm = T))
      MapY = c(MapY, mean(st_coordinates(M1[i,])[,'Y'], na.rm = T))
      m = which(tmp_map_place_names[i] == tmp_Place_names)
      Map2PlaceName = c(Map2PlaceName, m)
      I = c(I, i)
    }
  }

  M = M1[I,]

  save("M", "Map2PlaceName", "MapX", "MapY", file = "Scottish_Map.RData")
}



