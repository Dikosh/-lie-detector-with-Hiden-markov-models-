Tm<-10000 # Max time
M<-2 # Number of states
S<-seq(from=1, to=M) # States [1,2]
I<-rep(1./M, each=M) # Initial distribution [0.5, 0.5]
T<-matrix(c(1-4./Tm, 4./Tm, 4./Tm,1-4./Tm),2) #[0.9995 0.0005; 0.0005 0.9995]
P<-c(function(o){dnorm(o)}, function(o){dnorm(o,sd=1.1)}) # Vector of density functions for each state (N(0,1), N(0,1.21))

E<-function(P,s,o){P[[s]](o)} # Calculate probability of observed value o for state s.

Es<-function(E,P,o) {
  # Same for all states
  probs<-c()
  for(s in S) {
    probs<-append(probs, E(P,s,o))
  }
  return(probs)
}

viterbi<-function(S,I,T,P,E,O,tm) {
  delta<-I*Es(E,P,O[1]) # initialization for t=1
  prev_states<-cbind(rep(0,length(S))) # zeros
  for(t in 2:Tm) {
    m<-matrix(rep(delta,length(S)),length(S))*T # delta(s')*T[ss'] forall s,s'
    md<-apply(m,2,max) # search for max delta(s) by s' for all s
    max_delta<-apply(m,2,which.max)
    prev_states<-cbind(prev_states,max_delta)  
    delta<-md*Es(E,P,O[t]) # prepare next step delta(s)
  }
  return(list(delta=delta,ps=prev_states)) # return delta_Tm(s) and paths
}

restoreStates<-function(d,prev_states) {
  idx<-which.max(d)
  s<-c(idx)
  sz<-length(prev_states[1,])
  for(i in sz:2) {
    idx<-prev_states[idx,i]
    s<-append(s,idx,after=0)
  }
  return (as.vector(s))
}
