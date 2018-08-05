Generate_Randon_Network<-function(N,p){
  network<-matrix(0,N,N)
  colnames(network)<-1:N
  rownames(network)<-1:N
  index<-1:N
  for(i in 1:N){
    index_temp<-index[index!=i]
    for(j in index_temp){
      if(runif(1,min=0,max=1)<=p){
        network[i,j]<-1
      }
    } 
  }
  return(network)
}

Generate_Small_World_Network<-function(N,p,K){
  network_sw<-matrix(0,N,N)         #这个K表示每个节点都与它左右的各K/2个节点
  colnames(network_sw)<-1:N         #相连所以这个K要求是一个偶数           
  rownames(network_sw)<-1:N         # K must be an even number         
  index<-1:N
  for(i in 1:N){
    for(j in 1:(K/2)){
      if((i-j)<=0){
        network_sw[i,N+i-j]<-1
      }else{
        network_sw[i,i-j]<-1
      }
      if((i+j)>N){
        network_sw[i,i+j-N]<-1
      }else{
        network_sw[i,i+j]<-1
      }
    }
  }
  ###print("the original regular network_small world:")
  ###print(network_sw)
  ###Above codes generate a regular network, then we randomly reconnect some linkages
  for(i in 1:N){                                  #through following codes:
    index_origin<-index[network_sw[i,]==1]
    num_erase<-0
    for(j in index_origin){
      if(runif(1,0,1)<=p){
        network_sw[i,j]<-0
        num_erase<-num_erase+1
      }
    }
    if(num_erase>0){
      index_temp<-index[(network_sw[i,]==0)&(index!=i)]
      index_re<-sample(index_temp,num_erase)
      network_sw[i,index_re]<-1
    }
  }
  return(network_sw)
}  

Generate_Scale_Free_Network<-function(N,m0,m){ #note that m<m0. And since a lot of 
  network_sf<-matrix(0,N,N) #researches have already proved that the parameter gamma
  colnames(network_sf)<-1:N #network has nothing to do with the m0 and m of a scale 
  rownames(network_sf)<-1:N #free here in the Preferential Attachment mechanism and
  index<-1:N                #scale free networks generated in this way will always 
  for(i in 1:m0){           #have a gamma equaling to 3. So in this study we will
    index_temp<-index[(index!=i)&(index<=m0)] #gamma constantly equal to 3.If you want  
    for(j in index_temp){  #to only use scale free networks with generate some other 
      network_sf[i,j]<-1   #scale free network with different gamma, you need to adopt
    }                      #improved mechanisms and rewrite the the code accordingly.
  }
  ###print("the original core network_scalefree:")
  ###print(network_sf)
  for(i in (m0+1):N){
    out_dist<-apply(network_sf,1,sum)[1:i-1]
    in_dist<-apply(network_sf,2,sum)[1:i-1]                        
    out_temp<-out_dist/sum(out_dist)      #create PDF
    in_temp<-in_dist/sum(in_dist)
    for(k in 1:length(in_temp)){          #transfer the PDF to CDF
      out_dist[k]<-sum(out_temp[1:k])
      in_dist[k]<-sum(in_temp[1:k])
    }
    in_new<-vector("numeric",m)   #claim in_new and out_new.A very important thing 
    out_new<-vector("numeric",m)  #is to remember that in_new is linked to out_dist 
                                  #and out_new is linked to in_dist                
    for(q in 1:m){                       #realize in_new and out_new     
      icon_in<-0
      while(icon_in==0){
        in_ran<-runif(1,0,1)
        for(j in 1:length(out_dist)){
          if(in_ran<=out_dist[j]){
            if(sum(in_new==j)!=0){
              icon_in<-0
              break
            }else{
              in_new[q]<-j
              icon_in<-1
              break
            }
          }
        }
      }
    }             
    for(q in 1:m){                               
      icon_out<-0
      while(icon_out==0){
        out_ran<-runif(1,0,1)
        for(j in 1:length(in_dist)){
          if(out_ran<=in_dist[j]){
            if(sum(out_new==j)!=0){
              icon_out<-0
              break
            }else{
              out_new[q]<-j
              icon_out<-1
              break
            }
          }
        }
      }
    }                            #now we have newly updated in_new and out_new
    network_sf[in_new,i]<-1
    network_sf[i,out_new]<-1
  }
  return(network_sf)  
}


#one thing to be clarified: we consider the row name(on the left side of the matirx)
#as the starting point of a link and the coloum name(on the top of the matrix) to
#be the ending point of a link. So element (1,3) represents the link (1->3) while 
#element (3,1) represents the link (3->1).In conclusion, the whole structure
#information is now all stored in the matrix named "network"

#Then we store the whole personal information in another matrix named "info",
#actually we can also store these information in a list.


Generate_info<-function(network,A_M=0.8,A_IB=0.2,C=0.06){ #c means capital 
  N<-length(network[1,])
  info<-matrix(0,nrow=N,ncol=7)
  rownames(info)<-1:N
  colnames(info)<-c("In","Out","A_M","A_IB","L_IB","D","C")
  info[,"A_M"]<-A_M
  info[,"A_IB"]<-A_IB
  info[,"C"]<-C
  for(i in 1:N){
    info[i,"Out"]<-sum(network[i,])
    info[i,"In"]<-sum(network[,i])
  }
  Link_weight<<-network        #其实这个函数会产生两个东西：一个info，一个Link_Weight.
  Link_weight[,]<<-0           #这两个两个变量都很重要，但是我的设置是：让info作为返回
  for(i in 1:N){               #值从这个函数肛门处的return处返回出去，而让Link_Weight则
    if(info[i,"In"]!=0){       #从函数的腹部直接破肚而出，产生到工作台上去
      Link_weight[,i]<<-network[,i]*A_IB/(info[i,"In"]) 
    }else{
      info[i,"A_M"]<-1
      info[i,"A_IB"]<-0
    }
  }
  Link_weight<<-round(Link_weight,3)
  for(i in 1:N){
    info[i,"L_IB"]<-sum(Link_weight[i,])
    info[i,"D"]<-1-info[i,"C"]-info[i,"L_IB"]
  }
  return(round(info,3))
}

#到目前为止，我们已经有了产生一个随机网络具体样本的完备工具：也就是两个函数，
#分别用来产生network和info
#for now we have complete tools to generate a random network:two kinds of 
#functions, which are used to generate "network" and "info" respectively.
#-----------------------------------------------------------------------------#

Check_Bankrupt<-function(InfMatr,nodes,q=0.99){
  N<-length(InfMatr[,1])
  index<-1:N
  Bankrupt_Spectrum<<-Bankrupt_Spectrum*0     #每次使用试管之前先清洗上次的残余
  for(i in nodes){
    K<-InfMatr[i,"A_IB"]+q*InfMatr[i,"A_M"]-InfMatr[i,"D"]-InfMatr[i,"L_IB"]
    if(K<0){ 
      Bankrupt_Spectrum[i]<<-1
    }
  }
  Who_Bkrpt_This_Term<<-index[Bankrupt_Spectrum!=0]
}
#这个函数就不用设置返回值了，因为我把它的功能设置为处理一个工作台（函数外部）上
#的变量Who_Bkrpt_This_Term了，所以这个变量一直在工作台上，不用从函数里面吐出来
#或者说，这个函数会“返回（不是真正意义上从return处的返回哦）”两个变量：
#Bankrupt_Spectrum和Who_Bkrpt_This_Term，这两个变量无非是以不同的形式记录了同样
#的信息

Track_Exposure<-function(network,nodes){
  N<-length(network[1,])
  index<-1:N
  Exposure_Spectrum<<-Exposure_Spectrum*0     #使用前清洗试管
  for(i in as.character(nodes)){
    Exposure_Spectrum<<-Exposure_Spectrum+(network[i,]==1)
  }
  Who_Is_Exposed<<-index[Exposure_Spectrum!=0]
}
#这个函数也是“返回（不是真正意义上从return处的返回哦）”了两个变量：
#Exposure_Spectrum和Who_Is_Exposed，这两个变量同样是以不同的形式记录
#了相同的信息

#其实上述两个函数中，除了使用了Bankrupt_Spectrum和Exposure_Spectrum这两个试
#管容器以外，还使用了who_Is_Exposed和Who_Bkrpt_This_Term这两个容器，但是这两
#个容器是用后不用清洗的，由于他们是一次性乘装，所以他们每次盛装新内容时会自动
#覆盖掉前一次的残余，而Bankrupt_Spectrum和Exposure_Spectrum这两个容器则不然，
#由于这两个容器的填装是逐个元素地进行，所以如果不清洗，则上一次的残余会污染新
#的内容，另外要注意Bankrupt_Spectrum的每一项，要不是0，要不是1，分别表示没破产
#或破产了，而Exposure_Spectrum的每一项的值并不局限于0或1，而有可能是2、3、4...
#等等。这个数值的意义是：表示这一项的标号所对应的那个node在这一轮的筛查中被
#expose了几次。不过把spectrum中的信息转化为Who_Is_Exposed中的信息时，反正只要
#是不为零的就全部记为被expose了，不管它被expose了几次。

#以下是两个最重要的函数，分别用于更新info和network，其实有两种可行的对info或
#者network的更新方法：我原先想到的就是将每次破产的nodes彻底从两个矩阵里抹去，
#也就是真的把矩阵不断变小。但是后来想想这样真的去动矩阵的骨骼大小太麻烦了，要
#“将破产的点从矩阵上抹去”，其实只要把关于这个点对应的横行和竖行的记录全都清零
#即可,而且最让人满意的性质是：这样的清零，在图像实质上就是把这个node孤立出来，
#这样的话，在任何方面的效果上都完全等同于将这个node从图像上彻底抹去，因为它被
#孤立了，所以之后再用Check_Exposure去筛查的时候，也查不到它的头上。跟这个点不
#存在了效果是一样的。

#另外，在info阵中，一条横行的数据记录被清零，肯定是明显怪异的，可以直接识别为是
#既往的破产点，但是在network阵中，一个横行和竖行均全为0的点，可能一开始生成的时
#候就是一个孤立点，而不一定是事后破产所致。所以为了区分先天的孤立点和后天的破产
#点，我专门开列了一个破产名录Bankrupcy_List，用于记录后天陆续破产的点。


Update_network_info_list<-function(Who_Bkrpt_This_Term,network,InfMatr){
  N<-length(network[1,])
  index<-1:N
  for(i in Who_Bkrpt_This_Term){
    Out<-index[network[i,]==1]
    for(j in Out){
      InfMatr[j,"In"]<-InfMatr[j,"In"]-1
      InfMatr[j,"A_IB"]<-InfMatr[j,"A_IB"]-Link_weight[i,j]
      network[i,j]<-0
      Link_weight[i,j]<<-0
    }
    In<-index[network[,i]==1]
    for(k in In){
      InfMatr[k,"Out"]<-InfMatr[k,"Out"]-1
###   InfMatr[k,"L_IB"]<-InfMatr[k,"L_IB"]-Link_weight[k,i] //关于这行语句，很有
#讲究：
#    原先我自己的理解是：在这个简易的模型里，我们认为当一个点破产以后，这个
#点欠别别的点的钱就全部一分不还了，即recovery rate=0（破产清算的时候一分钱也不
#给债主）如果这样的话，那么对应的，这个点破产之时，欠了这个node钱的银行也不必再
#还它钱了，所以当一个node破产以后，所有它指向的nodes们，都要在各自的A_IB上扣除一
#笔钱（即借给这个破产的node的钱随着这个node的破产完全打水漂了，连个recovery也
#没有），而与此同时，指向这个破产node的nodes则赚了便宜，可以在它们各自的L_IB上
#扣除一笔钱即因为债主死了，所以欠它的钱也不用还了。当然，以上的操作是不切合实
#际的，是一种武断的简化操作，实际中一个node破产之前应该从所有欠了它钱的nodes那
#儿把那些钱都收回来，当然，另一方面，这个node破产以后，也不可能一点都不偿付它
#欠别的nodes的钱（即实际的recovery rate不可能是0），我删除这行语句，（以及添加
#其他地方的另外一些语句）就表示采用了后面一种更为现实的设计。不过在本代码中，我
#只是单方面删了这句语句而并未在其他地方添加相适应的另外语句，这是因为教授决定采
#用这样“跛脚的逻辑”以图与GK论文上的模型吻合。我觉得这种想法极其愚蠢（教授自己也
#承认了这样做本质上是不对的，如果真的要正确的话，就是应该按照我想的那样做，不过
#他就是为了迎合GK那篇原始论文上的做法，真是智障。）
      network[k,i]<-0
      Link_weight[k,i]<<-0
    }
    Bankrupcy_List<<-c(Bankrupcy_List,i) #这里是在更新破产名录Bankrupcy_List
    InfMatr[i,]<-0
  }
  return(list(network,InfMatr))  #这里也有两个返回值，但是我这不能在像之前在
}                                #Generate_info函数中对待info和Link_Weight那样
#将info用return返回，而把Link_Weight直接yong<<-剖腹产式地输出到工作台上。因为
#这里的network和InfMatr两者地位完全相当，是平权的。所以只能把他们俩攒成一个list
#然后整体返回这个list了，其实在Generate_info中，我也可以把info和Link_Weight攒
#成一个list来返回的。只是我当时没有这么做。

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#以下为主函数：Below is the main function:
#in the following codes, we use "NWK" to represent "network", so 
#"R_NWK"=random network; "SW_NWK"=small world network; "SF_NWK"=scale free network
#and we use "R_INF" to represent "info matrix of random network";
#"SW_INF" represents "info matrix of small world network"; "SF_INF" represents
#"info matrix of scale free network"
#-----------------------------------------------------------------------------
############################## FOR RANDOM NETWORK ############################
Q<-0.99
Capital<-0.06
a_m<-0.8
a_ib<-0.2
size<-35          ######可调(tunable)######                the size of the network        
Bankrupt_Spectrum<-vector("numeric",length=size)#这是两根公共试管，谁都可以来借用，
names(Bankrupt_Spectrum)<-1:size                #当然，要注意，每个来借这是两根公共
Exposure_Spectrum<-vector("numeric",length=size)#试管的人，在使用之前必须自己先清洗，
names(Exposure_Spectrum)<-1:size                #以防止上一次使用后的残留污染这次实验
#Bankrupt_Spectrum
#Exposure_Spectrum 
Bankrupcy_List<-vector("numeric",length=0)      #这就是上文提到的破产名录
#Bankrupcy_List
iteration_times<-100                             ######可调(tunable)######
Bankrupcy_rate_list<-vector("numeric",length=0)
#Bankrupcy_rate_list
AVE_DEGREE_LIST<-seq(0,9,0.1)                    ######可调(tunable)######
conditional_extent_of_contagion<-vector("numeric",length=0)
frequency_of_contagion<-vector("numeric",length=0)
contagion_threshold<-0.2                         ######可调(tunable)######

for(Ave_degree in AVE_DEGREE_LIST){
  Bankrupcy_rate_list<-vector("numeric",length=0)
  #Bankrupcy_rate_list
  for(i in 1:iteration_times){
    Bankrupt_Spectrum<-vector("numeric",length=size)
    names(Bankrupt_Spectrum)<-1:size
    #Bankrupt_Spectrum
    Exposure_Spectrum<-vector("numeric",length=size)
    names(Exposure_Spectrum)<-1:size
    #Exposure_Spectrum
    Bankrupcy_List<-vector("numeric",length=0)
    #Bankrupcy_List
    R_NWK<-Generate_Randon_Network(N=size,p=Ave_degree/(size-1))
    #R_NWK
    R_INF<-Generate_info(R_NWK,A_M=a_m,A_IB=a_ib,C=Capital)
    #Link_weight                      
    #R_INF                              
    initial_node<-sample(1:size,1)       
    #initial_node
    R_INF[initial_node,"A_M"]<-0
    #R_INF
    Check_Bankrupt(R_INF,initial_node,q=Q) 
    #Bankrupt_Spectrum             
    #Who_Bkrpt_This_Term           
    while(sum(Bankrupt_Spectrum)!=0||(length(Who_Bkrpt_This_Term)!=0)){
      Track_Exposure(R_NWK,Who_Bkrpt_This_Term)
      #print("Exposure Spectrum:")
      #print(Exposure_Spectrum)
      #print("Who is exposed:")
      #print(Who_Is_Exposed)
      result<-Update_network_info_list(Who_Bkrpt_This_Term,R_NWK,R_INF)
      R_NWK<-result[[1]]                                            
      R_INF<-result[[2]]
      #print("network:")  
      #print(R_NWK)
      #print("info:")
      #print(R_INF)
      Check_Bankrupt(R_INF,Who_Is_Exposed,q=Q)
      #print("bankrupt spectrum:")
      #print(Bankrupt_Spectrum)
      #print("who is bankrupt this term:")
      #print(Who_Bkrpt_This_Term)
    } 
    #print("bankrupt list:")
    #Bankrupcy_List
    #print("bankrupcy rate:")
    rate<-length(Bankrupcy_List)/size
    #rate
    Bankrupcy_rate_list<-c(Bankrupcy_rate_list,rate)
    #Bankrupcy_rate_list
  }
  pro<-mean(Bankrupcy_rate_list>contagion_threshold)
  frequency_of_contagion<-c(frequency_of_contagion,pro)
  exten<-mean(Bankrupcy_rate_list[Bankrupcy_rate_list>contagion_threshold])
  conditional_extent_of_contagion<-c(conditional_extent_of_contagion,exten)
}
plot(AVE_DEGREE_LIST,conditional_extent_of_contagion,xlab="average degree",
     ylab="probability",main="RANDOM NETWORK",pch=20,ylim=c(0,1),col="red")
lines(AVE_DEGREE_LIST,frequency_of_contagion,type="p",pch=19,cex=0.6,col="black")
legend("right",legend=c("Contagion frequency","Contagion extent(cdl)"),
       cex=0.6,fill=1:2)

#----------------------------------------------------------------------------------
########################## FOR SMALL WORLD NETWORK ################################
PPP<-0.4
Q<-0.99
Capital<-0.06
a_m<-0.8
a_ib<-0.2
size<-500          ######可调(tunable)######                the size of the network        
Bankrupt_Spectrum<-vector("numeric",length=size)#这是两根公共试管，谁都可以来借用，
names(Bankrupt_Spectrum)<-1:size                #当然，要注意，每个来借这是两根公共
Exposure_Spectrum<-vector("numeric",length=size)#试管的人，在使用之前必须自己先清洗，
names(Exposure_Spectrum)<-1:size                #以防止上一次使用后的残留污染这次实验
#Bankrupt_Spectrum
#Exposure_Spectrum 
Bankrupcy_List<-vector("numeric",length=0)      #这就是上文提到的破产名录
#Bankrupcy_List
iteration_times<-50                             ######可调(tunable)######
Bankrupcy_rate_list<-vector("numeric",length=0)
#Bankrupcy_rate_list
AVE_DEGREE_LIST<-seq(0,16,2)                    ######可调(tunable)######
conditional_extent_of_contagion<-vector("numeric",length=0)
frequency_of_contagion<-vector("numeric",length=0)
contagion_threshold<-0.002                    ######可调(tunable)######
#*******************************************************************************************
#*Try this: maintain all the other parameters unchanged,then change the contagion_threshold*
#*to be less than 0.002,say 0.001,very interestring phenomenon will be exibited,which is:  *
#*The features and attributes of the two curves will exchange!                             *
#*******************************************************************************************
for(Ave_degree in AVE_DEGREE_LIST){
  Bankrupcy_rate_list<-vector("numeric",length=0)
  #Bankrupcy_rate_list
  for(i in 1:iteration_times){
    Bankrupt_Spectrum<-vector("numeric",length=size)
    names(Bankrupt_Spectrum)<-1:size
    #Bankrupt_Spectrum
    Exposure_Spectrum<-vector("numeric",length=size)
    names(Exposure_Spectrum)<-1:size
    #Exposure_Spectrum
    Bankrupcy_List<-vector("numeric",length=0)  
    #Bankrupcy_List
    SW_NWK<-Generate_Small_World_Network(N=size,p=PPP,K=Ave_degree)
    #SW_NWK
    SW_INF<-Generate_info(SW_NWK,A_M=a_m,A_IB=a_ib,C=Capital)
    #Link_weight                      
    #SW_INF                              
    initial_node<-sample(1:size,1)       
    #initial_node
    SW_INF[initial_node,"A_M"]<-0
    #SW_INF
    Check_Bankrupt(SW_INF,initial_node,q=Q) 
    #Bankrupt_Spectrum             
    #Who_Bkrpt_This_Term           
    while(sum(Bankrupt_Spectrum)!=0||(length(Who_Bkrpt_This_Term)!=0)){
      Track_Exposure(SW_NWK,Who_Bkrpt_This_Term)
      #print("Exposure Spectrum:")
      #print(Exposure_Spectrum)
      #print("Who is exposed:")
      #print(Who_Is_Exposed)
      result<-Update_network_info_list(Who_Bkrpt_This_Term,SW_NWK,SW_INF)
      SW_NWK<-result[[1]]                                            
      SW_INF<-result[[2]]
      #print("network:")
      #print(SW_NWK)
      #print("info:")
      #print(SW_INF)
      Check_Bankrupt(SW_INF,Who_Is_Exposed,q=Q)
      #print("bankrupt spectrum:")
      #print(Bankrupt_Spectrum)
      #print("who is bankrupt this term:")
      #print(Who_Bkrpt_This_Term)
    }
    #print("bankrupt list:")
    #Bankrupcy_List
    #print("bankrupcy rate:")
    rate<-length(Bankrupcy_List)/size
    #rate
    Bankrupcy_rate_list<-c(Bankrupcy_rate_list,rate)
    #Bankrupcy_rate_list
  }
  pro<-mean(Bankrupcy_rate_list>contagion_threshold)
  frequency_of_contagion<-c(frequency_of_contagion,pro)
  exten<-mean(Bankrupcy_rate_list[Bankrupcy_rate_list>contagion_threshold])
  conditional_extent_of_contagion<-c(conditional_extent_of_contagion,exten)
}
plot(AVE_DEGREE_LIST,conditional_extent_of_contagion,xlab="average degree",
     ylab="probability",main="SMALL WORLD NETWORK",pch=20,ylim=c(0,1),col="red")
lines(AVE_DEGREE_LIST,frequency_of_contagion,type="p",pch=19,cex=0.6,col="black")
legend("right",legend=c("Contagion frequency","Contagion extent(cdl)"),
       cex=0.6,fill=1:2)

#Here we can observe a very interesting thing!!! The small world network seems to 
#exibit a very regular feature. Well, of course, the demonstration of such a rare
#and regular feature must also has something to do with the small world network
#generation mechanism, since we know that there are quite many mechanisms to generate
#a small world network.

#----------------------------------------------------------------------------------
############################# FOR SCALE FREE NETWORK #############################
M0<-10
Q<-0.99
Capital<-0.06
a_m<-0.8
a_ib<-0.2
size<-50          ######可调(tunable)######                the size of the network        
Bankrupt_Spectrum<-vector("numeric",length=size)#这是两根公共试管，谁都可以来借用，
names(Bankrupt_Spectrum)<-1:size                #当然，要注意，每个来借这是两根公共
Exposure_Spectrum<-vector("numeric",length=size)#试管的人，在使用之前必须自己先清洗，
names(Exposure_Spectrum)<-1:size                #以防止上一次使用后的残留污染这次实验
#Bankrupt_Spectrum
#Exposure_Spectrum 
Bankrupcy_List<-vector("numeric",length=0)      #这就是上文提到的破产名录
#Bankrupcy_List
iteration_times<-50                               ######可调(tunable)######
Bankrupcy_rate_list<-vector("numeric",length=0)
#Bankrupcy_rate_list
AVE_DEGREE_LIST<-seq(0,9,1)                       ######可调(tunable)######
conditional_extent_of_contagion<-vector("numeric",length=0)
frequency_of_contagion<-vector("numeric",length=0)
contagion_threshold<-0.03                     ######可调(tunable)######
           #****************************************************************
           #*The attributes of scale free network is very much analogous to*   
           #*those of small world network!                                 *
           #****************************************************************
for(Ave_degree in AVE_DEGREE_LIST){
  Bankrupcy_rate_list<-vector("numeric",length=0)
  #Bankrupcy_rate_list
  for(i in 1:iteration_times){
    Bankrupt_Spectrum<-vector("numeric",length=size)
    names(Bankrupt_Spectrum)<-1:size
    #Bankrupt_Spectrum
    Exposure_Spectrum<-vector("numeric",length=size)
    names(Exposure_Spectrum)<-1:size
    #Exposure_Spectrum
    Bankrupcy_List<-vector("numeric",length=0)  
    #Bankrupcy_List
    SF_NWK<-Generate_Scale_Free_Network(N=size,m0=M0,m=Ave_degree)
    #SF_NWK
    SF_INF<-Generate_info(SF_NWK,A_M=a_m,A_IB=a_ib,C=Capital)
    #Link_weight                      
    #SF_INF                              
    initial_node<-sample(1:size,1)       
    #initial_node
    SF_INF[initial_node,"A_M"]<-0
    #SF_INF
    Check_Bankrupt(SF_INF,initial_node,q=Q) 
    #Bankrupt_Spectrum             
    #Who_Bkrpt_This_Term           
    while(sum(Bankrupt_Spectrum)!=0||(length(Who_Bkrpt_This_Term)!=0)){
      Track_Exposure(SF_NWK,Who_Bkrpt_This_Term)
      #print("Exposure Spectrum:")
      #print(Exposure_Spectrum)
      #print("Who is exposed:")
      #print(Who_Is_Exposed)
      result<-Update_network_info_list(Who_Bkrpt_This_Term,SF_NWK,SF_INF)
      SF_NWK<-result[[1]]                                            
      SF_INF<-result[[2]]
      #print("network:")
      #print(SF_NWK)
      #print("info:")
      #print(SF_INF)
      Check_Bankrupt(SF_INF,Who_Is_Exposed,q=Q)
      #print("bankrupt spectrum:")
      #print(Bankrupt_Spectrum)
      #print("who is bankrupt this term:")
      #print(Who_Bkrpt_This_Term)
    }
    #print("bankrupt list:")
    #Bankrupcy_List
    #print("bankrupcy rate:")
    rate<-length(Bankrupcy_List)/size
    #rate
    Bankrupcy_rate_list<-c(Bankrupcy_rate_list,rate)
    #Bankrupcy_rate_list
  }
  pro<-mean(Bankrupcy_rate_list>contagion_threshold)
  frequency_of_contagion<-c(frequency_of_contagion,pro)
  exten<-mean(Bankrupcy_rate_list[Bankrupcy_rate_list>contagion_threshold])
  conditional_extent_of_contagion<-c(conditional_extent_of_contagion,exten)
}
plot(AVE_DEGREE_LIST,conditional_extent_of_contagion,xlab="average degree",
     ylab="probability",main="SCALE FREE NETWORK",pch=20,ylim=c(0,1),col="red")
lines(AVE_DEGREE_LIST,frequency_of_contagion,type="p",pch=19,cex=0.6,col="black")
legend("right",legend=c("Contagion frequency","Contagion extent(cdl)"),
       cex=0.6,fill=1:2)

