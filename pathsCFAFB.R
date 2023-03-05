library(natverse)
library(neuprintr); library(hemibrainr); library(nat); library(natverse); library(nat.jrcbrains); library(fafbseg);
library(reticulate)
library(R.matlab)
# library(tidyverse)
library(pheatmap)
library(igraph)
library(pracma)
library(visNetwork)
library(zoo)


### METHODS
getDownstreamsNRFAFB <- function(ids=c(), nLayers=2, connStrength=c(5,10)) {
  s_ConnT=c(list())
  s_Conn=c()
  sConn=list()
  counter=c()
  itr=2
  
   for (i in 1:nLayers) { # Loop through Number of desired layers
    conn_N=flywire_partner_summary(rootids=ids,partners = 'outputs',threshold = connStrength[i])
    
    nconn_N=conn_N['post_id'][[1]]

    if (i==1) {
      
      ids=unique(nconn_N)# Repeat process if i ~= nLayers
      counterT=rep(i,length(ids))
      counter=c(counter,counterT)
      s_Conn=c(s_Conn,ids)
      sConn[[itr]]=conn_N
      itr=itr+1
    } else {
      
      P=setdiff(unique(nconn_N),s_Conn) # Extract connections unique to this layer with respect to previous layers
      counterT=rep(i,length(P))
      counter=c(counter,counterT)
      s_Conn=c(s_Conn,P)
      sConn[[itr]]=conn_N
      
      ids=P# Repeat process if i ~= nLayers
      itr=itr+1
    }
   }
  s_Conn=data.frame(s_Conn,counter)
  sConn[[1]]=s_Conn
  return(sConn)
}

# clN=read.csv("C:\\Users\\LabAdmin\\Downloads\\search_results_cl062.csv")$root_id
# cLN=read.csv("C:\\Users\\LabAdmin\\Downloads\\root_ids_cl062_.csv")
# clN=c('720575940660182145','720575940626885124','720575940634371941','720575940608082990','720575940630037073','720575940637079667','720575940628333788','720575940606687957','720575940614437462','720575940620424472','720575940626869403','720575940632011420','720575940620123613','720575940637319134','720575940602120928','720575940641972384','720575940613034602','720575940604445868','720575940619688688','720575940617362772','720575940637175349','720575940617034398')

### *EDIT HERE*

### INPUT NEUROPIL
clNR=c('720575940626885124','720575940634371941','720575940608082990','720575940630037073'
       ,'720575940637079667','720575940606687957','720575940614437462','720575940641972384',
       '720575940604445868','720575940617034398')
clNL=c('720575940628333788','720575940620424472','720575940626869403','720575940632011420',
       '720575940620123613','720575940637319134','720575940619688688','720575940637175349')
clN=c(clNR,clNL)

### TARGET NEUROPIL
DNfw=strsplit(readLines("C:\\Users\\LabAdmin\\Downloads\\root_ids_class_equal_descending.txt"),',')[[1]]

### VECTOR OF THRESHOLDS PER LAYER
Thr=c(10,15)
### NUMBER OF LAYERS TO TRAVERSE 
nLayers=2

# inpN=c()
# for (i in clN){
#   inpN=c(inpN,toString(i))
# }

### SET INPUT NEUROPIL NEURONS
inpN=clN
# inpN=clNR
# inpN=clNL


### ASSEMBLE ADJACENCY MATRIX BY FAFBSEG'S FUNCTION OR DOWNSTREAM CONNECTIONS
byConn=TRUE

### PATH LENGTHS TO CONSIDER (LEAVE EMPTY TO CONSIDER ALL)
pL=c(1:2)

### MAIN

sConnAll=getDownstreamsNRFAFB(inpN,nLayers,Thr)
sConnA=sConnAll[[1]]$s_Conn

if (byConn==TRUE) {

adjFW=matrix(0,length(sConnA),length(sConnA))
rownames(adjFW)=sConnA
colnames(adjFW)=sConnA


for (k in 2:length(sConnAll)) {
    for (d in 1:nrow(sConnAll[[k]])) {
      adjFW[which(rownames(adjFW)%in%sConnAll[[k]]$query[d]),which(colnames(adjFW)%in%sConnAll[[k]]$post_id[d])]=sConnAll[[k]]$weight[d]
  }
}


} else {

adjFW=flywire_adjacency_matrix(c(inpN,sConnA))
}

G=graph_from_adjacency_matrix(adjFW)



fromt=c()
from=list()
tot=c()
to=list()
wtt=c()
wt=list()
wtc=c()
dpt=list()
outN=rownames(adjFW)[which(rownames(adjFW)%in%DNfw)]
itr=1

if (length(pL)==0){
for (k in which(rownames(adjFW)%in%inpN)){
  allP=shortest_paths(G,k,to=which(rownames(adjFW)%in%DNfw))
  allP=list(allP$vpath)
for (i in allP[[1]]) {
  fromt=c()
  tot=c()
  wtt=c()
  dptt=c()
  if (rownames(adjFW)[i[length(i)]]%in%outN) {
      for (n in 1:(length(i)-1)){
      fromt=c(fromt,rownames(adjFW)[i[n]])
      tot=c(tot,colnames(adjFW)[i[n+1]])
      wtt=c(wtt,1/adjFW[i[n],i[n+1]])
      }
 
  dptt=1:(length(i)-1)
  dpt[[itr]]=dptt
   wtc=c(wtc,(sum(wtt)))
   from[[itr]]=fromt
   to[[itr]]=tot
   wt[[itr]]=wtt
   itr=itr+1
  } 
  }
}
} else {
  for (k in which(rownames(adjFW)%in%inpN)){
    allP=shortest_paths(G,k,to=which(rownames(adjFW)%in%DNfw))
    allP=list(allP$vpath)
    for (i in allP[[1]]) {
      fromt=c()
      tot=c()
      wtt=c()
      dptt=c()
      if (rownames(adjFW)[i[length(i)]]%in%outN) {
       if ((length(i)-1)%in%pL) { 
         for (n in 1:(length(i)-1)){
          
          fromt=c(fromt,rownames(adjFW)[i[n]])
          tot=c(tot,colnames(adjFW)[i[n+1]])
          wtt=c(wtt,1/adjFW[i[n],i[n+1]])
          }
       
      dptt=1:(length(i)-1)
      dpt[[itr]]=dptt
      wtc=c(wtc,(sum(wtt)))
      from[[itr]]=fromt
      to[[itr]]=tot
      wt[[itr]]=wtt
      itr=itr+1
        } 
      }
    }
  } 
}

pathsn=data.frame()
fromn=c()
ton=c()
wtn=c()
dptn=c()
for (i in order(wtc)) {
  fromn=c(fromn,from[[i]])
  ton=c(ton,to[[i]])
  wtn=c(wtn,1/wt[[i]])
  dptn=c(dptn,dpt[[i]])
  
}
pathsn=data.frame(fromn,ton,wtn,dptn)

distComp=data.frame();

colnames(pathsn)=c('from','to','weight','depth')

itr=1

allPathsO=list(pathsn)

for (j in allPathsO){
  distCompt=data.frame()
  for (i in 1:(nrow(j))) {
    if (i==nrow(j)) {
    distCompt[itr,1]=sum(1/(j$weight[(i-j$depth[i]+1):i]));
    distCompt[itr,2]=(j$to[i])
    distCompt[itr,3]=(j$from[i-j$depth[i]+1])
    }
    else if ((j$to[i]%in%outN)&(j$depth[i+1]==1)){
    distCompt[itr,1]=sum(1/(j$weight[(i-j$depth[i]+1):i]));
    distCompt[itr,2]=(j$to[i])
    distCompt[itr,3]=(j$from[i-j$depth[i]+1])
    itr=itr+1
    }
    
  }
  distComp=distCompt
}
effDF=data.frame(matrix(0,nrow=length(inpN),ncol=(length(unique((outN))))))
rownames(effDF)=c(unique((unlist(inpN))))
colnames(effDF)=unique((outN))



for (i in 1:nrow(distComp)) {
  effDF[which(rownames(effDF)==distComp[i,3]),which(colnames(effDF)==distComp[i,2])]=sum(effDF[which(rownames(effDF)==distComp[i,3]),which(colnames(effDF)==distComp[i,2])],1/distComp[i,1])
}

effDF=1/effDF

from=c()
to=c()
label=c()
for (i in 1:nrow(effDF)){
  for (j in 1:ncol(effDF)) {
    if (effDF[i,j]=='Inf') {
      j=j+1
    } else{
      from=c(from,rownames(effDF)[i])
      to=c(to,colnames(effDF)[j])
      label=c(label,effDF[i,j])
    }
  }
}

nodes=data.frame(id=c(unique(from),unique(to)),label=c(unique(from),unique(to)),level=c(rep(1,length(unique(from))),rep(2,length(unique(to)))))
edges=data.frame(from=from,to=to,value=1/label)

try(visNetwork(nodes, edges,width = '100%',height = ) %>% 
      visEdges(shadow = TRUE,
               arrows ='to',
               selectionWidth = 5,
               color = list(color = "lightblue", highlight = "red"))  %>% 
      visOptions(selectedBy = "label", 
                 highlightNearest = TRUE, 
                 nodesIdSelection = TRUE) %>%
      visHierarchicalLayout(levelSeparation = 500) %>%
      visPhysics(stabilization = FALSE,repulsion = FALSE))

write.csv(effDF,"C:\\Users\\LabAdmin\\Desktop\\Connectomics\\Deven Connectomics\\effDFallDNsLRFAFBpL12v2.csv", row.names = TRUE)

View(effDF)

RtoWrange<-colorRampPalette(c('white','red' ))
col=RtoWrange(100)
pheatmap(1/effDF,
         color = c(col),cluster_rows = FALSE,cluster_cols = FALSE,
         annotation_names_row = TRUE,
         annotation_names_col = TRUE,
         show_rownames = TRUE,show_colnames = TRUE,
         fontsize=5)


