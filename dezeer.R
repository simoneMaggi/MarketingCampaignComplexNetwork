library(igraph)
####### caricare dataset
data <- read.csv("deezer_europe/deezer_europe_edges.csv", sep = ',', header = TRUE)

data <- as.matrix(data)+1

G = graph_from_edgelist(data, directed = FALSE)

A <- as_adjacency_matrix(G)

n = vcount(G)
m = ecount(G)

# SMALL WORLD


global_trans <- transitivity(graph = G, type = "global")
local_trans <- transitivity(graph = G, type = "local")

rho = m/(n*(n-1)/2)

#avg shortest path
#2 minutes
distribution_geodetics <- distance_table(G)$res 
plot(log(distribution_geodetics))

#avg_path_lenght<- mean_distance(G, directed = FALSE)
avg_path_lenght <- sum(
  distribution_geodetics/sum(distribution_geodetics)*
    1:length(distribution_geodetics))

log(n)/log(2*m/n)

### check if it has hubs and the degree distribution follow a power law

degr_distr <- degree_distribution(graph = G)

lin_reg <- lm(log(degr_distr)~log(1:length(degr_distr)), subset = degr_distr>0)
summary(lin_reg)
plot(degr_distr, log = "xy")
abline(a = lin_reg$coefficients[1], b=lin_reg$coefficients[2], col = 'red')
C = exp(lin_reg$coefficients[1])
b = lin_reg$coefficients[2]

plot(degr_distr)
axe = 1:length(degr_distr)
lines(axe, C*(axe)^b, col='red')

#######  degree correlation

degree_correlation_func <- knn(G)$knnk

maxdegree <- length(degr_distr)
k_square <- sum(degr_distr * (1:maxdegree)^2)
k_mean <- sum(degr_distr * 1:maxdegree)
k_neutral <- k_square/k_mean

v <- 1:length(degree_correlation_func)

degr_cor_fit <- lm(degree_correlation_func ~ v, 
                   subset = is.nan(degree_correlation_func) == FALSE)

summary(degr_cor_fit)

plot(degree_correlation_func)
abline(h = k_neutral, col = "red")
abline(a = degr_cor_fit$coefficients[1], b=degr_cor_fit$coefficients[2], 
       col ="green")
assortativity_degree(G)

############   rich club effect

#we need to find the density in the subgraph
#induced by the fraction r of the most popular verteces


degrees <- Matrix::colSums(A)
#this vector contains the id nodes in decreasing order w.r.t. degree
decreas_nodes_degree<-order(degrees, decreasing = TRUE)

density_induced<- function(nodes_id){
  #computes the density of the subgraph induced by nodes_id
  n_d <- length(nodes_id)
  return(sum(A[nodes_id,nodes_id])/(n_d*(n_d-1)))
}

r <- seq(from=0, to=1, by=0.001)
phi <- sapply(r,function(x){
  nodes_id <- decreas_nodes_degree[1:round(n*x,0)]
  return(density_induced(nodes_id))
})

plot(phi, log="x")
#a clear rich club effect


# HOMOPHILY

# se fosse una rete random, senza omofilia(la rete random viene generata riposizionando
# casualmente gli edge tra i nodi), la probabilità di avere 
# un edge tra 1 e 2 è P = 2p(1-p) ovvero la probabilità di prendere
# un nodo dal gruppo 1 e uno dal gruppo 2 o viceversa ( da qui il 2)

# chiama Xij la r.v. associata ad ogni edge che è 1 se i-j è un link tra G1 e G2,
# e 0 altrimenti. Se la rete non avesse omofilia u0 = E[Xij] = 2p(1-p)
# TEST H0 : u = "media empirica di Xij" = u0 H1: u < u0
# call y = sqrt(n)(u-u0), as n -> inf, under H0 y is normal with mean = 0,
# var = u0(1-u0). Quindi trovo theta t. c. P(y < theta | H0) = alpha
# se y_data < theta allora posso accettare H1


sex <- read.csv("deezer_europe/sex.csv")
sex <- sex[,2]
table(sex)#male 0, woman 1


#number of links between groups
p <- sum(sex == 0)/n # fraction in group 1 
groups <- sex
vertices = 1:n
A = as_adjacency_matrix(G)

u0 = 2*p*(1-p)# fraction of edges in random network

#actual fraction of nodes
u=sum(A[vertices[groups==0],vertices[groups == 1]])/m
u<u0
# smaller but how much?
y = sqrt(m)*(u - u0)
#by CLT y is approx gaussian with mean 0 and variance u0(1-u0)
alpha = 0.05
pnorm(y, sd=sqrt(u0*(1-u0)))
theta = qnorm(alpha, mean = 0, sd = sqrt(u0*(1-u0)))
theta0 = theta/sqrt(m) + u0
u<theta0
# u minore della soglia thetha0 quindi possiamo accettare homophily



##### simulazione diffusione di opinione


degrees_inv <- 1/Matrix::colSums(A)

D_inv <- Diagonal(x = degrees_inv)
# deve essere una sparse matrix sennò pesa troppo

simulation <- function(x, q){  
  status <- rep(FALSE, n) # a vector that link nodes to conviced (0 or 1)
  if( length(x) == 1){
    S = sample(x = 1:n, size = x, replace = FALSE) # starting set of nodes
  }else{
    S = x
  }
  n_v<-length(S)
  status[S] <- rep(TRUE,n_v)
  n_convinced <- c(0, n_v)
  stp <- 2
  
  while (n_convinced[stp]-n_convinced[stp-1] > 0){
    #a node can be convinced if the fraction of neighbours convinced p is 
    # p*a>(1-p)*b => p > b/(a+b)=rho. 
    # status*A*D_inv gives for each node the value p
    print(paste("step ",stp," convinti ",sum(status)))
    
    status <- matrix(ifelse(status%*%A%*%D_inv > q, TRUE, FALSE), nrow = 1) | status 
    n_convinced <- c(n_convinced, sum(status))
    stp <- stp+1
  }
  
  return(n_convinced)
}

k=200
q=1/3


simulation(k,q)

# choosing as starting nodes the k most popular

simulation(decreas_nodes_degree[1:k],q)

#choosing the k most close 

close <- closeness(G)
decreas_nodes_close <- order(close, decreasing = TRUE)

simulation(decreas_nodes_close[1:k],q)

##### k-shell algorithm

kd <- coreness(G)

table(kd)
degrees[which(kd == 12)]

S = c(which(kd==12),which(kd == 11),which(kd == 10))
simulation(S,q)
