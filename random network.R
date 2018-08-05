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
  network_sw<-matrix(0,N,N)         #���K��ʾÿ���ڵ㶼�������ҵĸ�K/2���ڵ�
  colnames(network_sw)<-1:N         #�����������KҪ����һ��ż��           
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
  Link_weight<<-network        #��ʵ����������������������һ��info��һ��Link_Weight.
  Link_weight[,]<<-0           #��������������������Ҫ�������ҵ������ǣ���info��Ϊ����
  for(i in 1:N){               #ֵ������������Ŵ���return�����س�ȥ������Link_Weight��
    if(info[i,"In"]!=0){       #�Ӻ����ĸ���ֱ���ƶǶ���������������̨��ȥ
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

#��ĿǰΪֹ�������Ѿ����˲���һ�������������������걸���ߣ�Ҳ��������������
#�ֱ���������network��info
#for now we have complete tools to generate a random network:two kinds of 
#functions, which are used to generate "network" and "info" respectively.
#-----------------------------------------------------------------------------#

Check_Bankrupt<-function(InfMatr,nodes,q=0.99){
  N<-length(InfMatr[,1])
  index<-1:N
  Bankrupt_Spectrum<<-Bankrupt_Spectrum*0     #ÿ��ʹ���Թ�֮ǰ����ϴ�ϴεĲ���
  for(i in nodes){
    K<-InfMatr[i,"A_IB"]+q*InfMatr[i,"A_M"]-InfMatr[i,"D"]-InfMatr[i,"L_IB"]
    if(K<0){ 
      Bankrupt_Spectrum[i]<<-1
    }
  }
  Who_Bkrpt_This_Term<<-index[Bankrupt_Spectrum!=0]
}
#��������Ͳ������÷���ֵ�ˣ���Ϊ�Ұ����Ĺ�������Ϊ����һ������̨�������ⲿ����
#�ı���Who_Bkrpt_This_Term�ˣ������������һֱ�ڹ���̨�ϣ����ôӺ��������³���
#����˵����������ᡰ���أ��������������ϴ�return���ķ���Ŷ��������������
#Bankrupt_Spectrum��Who_Bkrpt_This_Term�������������޷����Բ�ͬ����ʽ��¼��ͬ��
#����Ϣ

Track_Exposure<-function(network,nodes){
  N<-length(network[1,])
  index<-1:N
  Exposure_Spectrum<<-Exposure_Spectrum*0     #ʹ��ǰ��ϴ�Թ�
  for(i in as.character(nodes)){
    Exposure_Spectrum<<-Exposure_Spectrum+(network[i,]==1)
  }
  Who_Is_Exposed<<-index[Exposure_Spectrum!=0]
}
#�������Ҳ�ǡ����أ��������������ϴ�return���ķ���Ŷ����������������
#Exposure_Spectrum��Who_Is_Exposed������������ͬ�����Բ�ͬ����ʽ��¼
#����ͬ����Ϣ

#��ʵ�������������У�����ʹ����Bankrupt_Spectrum��Exposure_Spectrum��������
#���������⣬��ʹ����who_Is_Exposed��Who_Bkrpt_This_Term��������������������
#���������ú�����ϴ�ģ�����������һ���Գ�װ����������ÿ��ʢװ������ʱ���Զ�
#���ǵ�ǰһ�εĲ��࣬��Bankrupt_Spectrum��Exposure_Spectrum������������Ȼ��
#������������������װ�����Ԫ�صؽ��У������������ϴ������һ�εĲ������Ⱦ��
#�����ݣ�����Ҫע��Bankrupt_Spectrum��ÿһ�Ҫ����0��Ҫ����1���ֱ��ʾû�Ʋ�
#���Ʋ��ˣ���Exposure_Spectrum��ÿһ���ֵ����������0��1�����п�����2��3��4...
#�ȵȡ������ֵ�������ǣ���ʾ��һ��ı������Ӧ���Ǹ�node����һ�ֵ�ɸ���б�
#expose�˼��Ρ�������spectrum�е���Ϣת��ΪWho_Is_Exposed�е���Ϣʱ������ֻҪ
#�ǲ�Ϊ��ľ�ȫ����Ϊ��expose�ˣ���������expose�˼��Ρ�

#��������������Ҫ�ĺ������ֱ����ڸ���info��network����ʵ�����ֿ��еĶ�info��
#��network�ĸ��·�������ԭ���뵽�ľ��ǽ�ÿ���Ʋ���nodes���״�����������Ĩȥ��
#Ҳ������İѾ��󲻶ϱ�С�����Ǻ��������������ȥ������Ĺ�����С̫�鷳�ˣ�Ҫ
#�����Ʋ��ĵ�Ӿ�����Ĩȥ������ʵֻҪ�ѹ���������Ӧ�ĺ��к����еļ�¼ȫ������
#����,��������������������ǣ����������㣬��ͼ��ʵ���Ͼ��ǰ����node����������
#�����Ļ������κη����Ч���϶���ȫ��ͬ�ڽ����node��ͼ���ϳ���Ĩȥ����Ϊ����
#�����ˣ�����֮������Check_Exposureȥɸ���ʱ��Ҳ�鲻������ͷ�ϡ�������㲻
#������Ч����һ���ġ�

#���⣬��info���У�һ�����е����ݼ�¼�����㣬�϶������Թ���ģ�����ֱ��ʶ��Ϊ��
#�������Ʋ��㣬������network���У�һ�����к����о�ȫΪ0�ĵ㣬����һ��ʼ���ɵ�ʱ
#�����һ�������㣬����һ�����º��Ʋ����¡�����Ϊ����������Ĺ�����ͺ�����Ʋ�
#�㣬��ר�ſ�����һ���Ʋ���¼Bankrupcy_List�����ڼ�¼����½���Ʋ��ĵ㡣


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
###   InfMatr[k,"L_IB"]<-InfMatr[k,"L_IB"]-Link_weight[k,i] //����������䣬����
#������
#    ԭ�����Լ��������ǣ���������׵�ģ���������Ϊ��һ�����Ʋ��Ժ����
#��Ƿ���ĵ��Ǯ��ȫ��һ�ֲ����ˣ���recovery rate=0���Ʋ������ʱ��һ��ǮҲ��
#��ծ������������Ļ�����ô��Ӧ�ģ�������Ʋ�֮ʱ��Ƿ�����nodeǮ������Ҳ������
#����Ǯ�ˣ����Ե�һ��node�Ʋ��Ժ�������ָ���nodes�ǣ���Ҫ�ڸ��Ե�A_IB�Ͽ۳�һ
#��Ǯ�����������Ʋ���node��Ǯ�������node���Ʋ���ȫ��ˮƯ�ˣ�����recoveryҲ
#û�У��������ͬʱ��ָ������Ʋ�node��nodes��׬�˱��ˣ����������Ǹ��Ե�L_IB��
#�۳�һ��Ǯ����Ϊծ�����ˣ�����Ƿ����ǮҲ���û��ˡ���Ȼ�����ϵĲ����ǲ��к�ʵ
#�ʵģ���һ����ϵļ򻯲�����ʵ����һ��node�Ʋ�֮ǰӦ�ô�����Ƿ����Ǯ��nodes��
#������ЩǮ���ջ�������Ȼ����һ���棬���node�Ʋ��Ժ�Ҳ������һ�㶼��������
#Ƿ���nodes��Ǯ����ʵ�ʵ�recovery rate��������0������ɾ��������䣬���Լ�����
#�����ط�������һЩ��䣩�ͱ�ʾ�����˺���һ�ָ�Ϊ��ʵ����ơ������ڱ������У���
#ֻ�ǵ�����ɾ�����������δ�������ط���������Ӧ��������䣬������Ϊ���ھ�����
#���������˽ŵ��߼�����ͼ��GK�����ϵ�ģ���Ǻϡ��Ҿ��������뷨�����޴��������Լ�Ҳ
#�������������������ǲ��Եģ�������Ҫ��ȷ�Ļ�������Ӧ�ð��������������������
#������Ϊ��ӭ��GK��ƪԭʼ�����ϵ��������������ϡ���
      network[k,i]<-0
      Link_weight[k,i]<<-0
    }
    Bankrupcy_List<<-c(Bankrupcy_List,i) #�������ڸ����Ʋ���¼Bankrupcy_List
    InfMatr[i,]<-0
  }
  return(list(network,InfMatr))  #����Ҳ����������ֵ���������ⲻ������֮ǰ��
}                                #Generate_info�����жԴ�info��Link_Weight����
#��info��return���أ�����Link_Weightֱ��yong<<-�ʸ���ʽ�����������̨�ϡ���Ϊ
#�����network��InfMatr���ߵ�λ��ȫ�൱����ƽȨ�ġ�����ֻ�ܰ��������ܳ�һ��list
#Ȼ�����巵�����list�ˣ���ʵ��Generate_info�У���Ҳ���԰�info��Link_Weight��
#��һ��list�����صġ�ֻ���ҵ�ʱû����ô����

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#����Ϊ��������Below is the main function:
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
size<-35          ######�ɵ�(tunable)######                the size of the network        
Bankrupt_Spectrum<-vector("numeric",length=size)#�������������Թܣ�˭�����������ã�
names(Bankrupt_Spectrum)<-1:size                #��Ȼ��Ҫע�⣬ÿ������������������
Exposure_Spectrum<-vector("numeric",length=size)#�Թܵ��ˣ���ʹ��֮ǰ�����Լ�����ϴ��
names(Exposure_Spectrum)<-1:size                #�Է�ֹ��һ��ʹ�ú�Ĳ�����Ⱦ���ʵ��
#Bankrupt_Spectrum
#Exposure_Spectrum 
Bankrupcy_List<-vector("numeric",length=0)      #����������ᵽ���Ʋ���¼
#Bankrupcy_List
iteration_times<-100                             ######�ɵ�(tunable)######
Bankrupcy_rate_list<-vector("numeric",length=0)
#Bankrupcy_rate_list
AVE_DEGREE_LIST<-seq(0,9,0.1)                    ######�ɵ�(tunable)######
conditional_extent_of_contagion<-vector("numeric",length=0)
frequency_of_contagion<-vector("numeric",length=0)
contagion_threshold<-0.2                         ######�ɵ�(tunable)######

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
size<-500          ######�ɵ�(tunable)######                the size of the network        
Bankrupt_Spectrum<-vector("numeric",length=size)#�������������Թܣ�˭�����������ã�
names(Bankrupt_Spectrum)<-1:size                #��Ȼ��Ҫע�⣬ÿ������������������
Exposure_Spectrum<-vector("numeric",length=size)#�Թܵ��ˣ���ʹ��֮ǰ�����Լ�����ϴ��
names(Exposure_Spectrum)<-1:size                #�Է�ֹ��һ��ʹ�ú�Ĳ�����Ⱦ���ʵ��
#Bankrupt_Spectrum
#Exposure_Spectrum 
Bankrupcy_List<-vector("numeric",length=0)      #����������ᵽ���Ʋ���¼
#Bankrupcy_List
iteration_times<-50                             ######�ɵ�(tunable)######
Bankrupcy_rate_list<-vector("numeric",length=0)
#Bankrupcy_rate_list
AVE_DEGREE_LIST<-seq(0,16,2)                    ######�ɵ�(tunable)######
conditional_extent_of_contagion<-vector("numeric",length=0)
frequency_of_contagion<-vector("numeric",length=0)
contagion_threshold<-0.002                    ######�ɵ�(tunable)######
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
size<-50          ######�ɵ�(tunable)######                the size of the network        
Bankrupt_Spectrum<-vector("numeric",length=size)#�������������Թܣ�˭�����������ã�
names(Bankrupt_Spectrum)<-1:size                #��Ȼ��Ҫע�⣬ÿ������������������
Exposure_Spectrum<-vector("numeric",length=size)#�Թܵ��ˣ���ʹ��֮ǰ�����Լ�����ϴ��
names(Exposure_Spectrum)<-1:size                #�Է�ֹ��һ��ʹ�ú�Ĳ�����Ⱦ���ʵ��
#Bankrupt_Spectrum
#Exposure_Spectrum 
Bankrupcy_List<-vector("numeric",length=0)      #����������ᵽ���Ʋ���¼
#Bankrupcy_List
iteration_times<-50                               ######�ɵ�(tunable)######
Bankrupcy_rate_list<-vector("numeric",length=0)
#Bankrupcy_rate_list
AVE_DEGREE_LIST<-seq(0,9,1)                       ######�ɵ�(tunable)######
conditional_extent_of_contagion<-vector("numeric",length=0)
frequency_of_contagion<-vector("numeric",length=0)
contagion_threshold<-0.03                     ######�ɵ�(tunable)######
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
