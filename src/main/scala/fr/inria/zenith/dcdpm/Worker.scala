package fr.inria.zenith.dcdpm

import scala.collection.mutable.ArrayBuffer
import org.apache.commons.math3.distribution.MultivariateNormalDistribution
import org.apache.commons.math3.distribution.BetaDistribution
import org.apache.commons.math3.distribution.GammaDistribution

class Worker(workerId : Int, val data : ArrayBuffer[Point], gamma : Double, val E : Double, val S : Double, var betaU : Double) {
  // dimension
  val dim = data(0).vector.length
  
  val nbPoint : Int = data.size
  // assignements
  var c : Map[Int, Int] = Map()
  
  var clusterVector = new ArrayBuffer[Cluster]
  
  val m = new Array[Double](dim)
    for(i <- 0 to dim - 1)
      m(i) = 0
  // Variance between centers
  var sigma1 = Array.ofDim[Double](dim, dim)
  for(i <- 0 to dim - 1)
    for(j <- 0 to dim - 1)
      if(i == j)
        sigma1(i)(j)= E
      else
        sigma1(i)(j) = 0
  // Variance in each cluster
  var sigma2 = Array.ofDim[Double](dim, dim)
  for(i <- 0 to dim - 1)
    for(j <- 0 to dim - 1)
      if(i == j)
        sigma2(i)(j)= S
      else
        sigma2(i)(j) = 0
        
  var alpha : Double = 1
  
  val H : Int = 3
        
  def getAlpha = alpha
  def getC = c
  
  def getCluster(id : Int) = {
    
    val index = clusterVector.indexWhere(cluster => (cluster.id == id))
    clusterVector(index)
    
  }
  
  def addCluster(index : Int, c : Cluster) = {
    
    clusterVector.insert(index, c)
    
  }
 
  def removePointFromCluster(clusterId : Int)={
    
    val index = clusterVector.indexWhere(cluster => (cluster.id == clusterId))
    
    clusterVector(index).size -= 1
    
  }
  
   def startWithClustering(globalClusters : ArrayBuffer[((Array[Double], Double), ArrayBuffer[(Int, Int)])], u : Double) = {
   
    betaU = u
     
    clusterVector.clear
    var newC : Map[Int, Int] = Map()
    
     var id = 0
    
    for (cluster <- globalClusters){
      
    val a = new Array[Double](dim)
    for(i <- 0 to dim - 1)
      a(i) = 0
      
      val clust = new Cluster(id, 0, a, cluster._1._1, cluster._1._2, false)
      
      
      for(subC <- cluster._2.filter(_._1 == workerId))
        for(ci <- c.filter({case(k,v) => v == subC._2})){
          newC += ci._1 -> id
          clust.addPoint()
        }
      
      clusterVector += clust
      id+=1
      
    }
    
    c = newC
        
    this
        
  }
  
  def initializeWithClustering(initial : ArrayBuffer[(Array[Double], Double)], u : Double) = {
    
    betaU = u
    clusterVector.clear
    c = Map()
    
    var id = 0
    
    for (tuple <- initial){

    val a = new Array[Double](dim)
    for(i <- 0 to dim - 1)
      a(i) = 0
      
      val clust = new Cluster(id, 0, a, tuple._1, tuple._2, false)
      clusterVector += clust
      id += 1
      
    }
    
    for(y <- data){
      
      val clusterId : Int = max(calculateProbaWithLog(y))
        
        c += (y.id -> clusterId)
        
        getCluster(clusterId).addPoint()
     
    }
        
    this

  }
  
  def gibbsSampling()={
      
    for(iter <- 0 to 9){
            
      // for each data point
      for(y <- data){
        // remove it from its cluster
        removePointFromCluster(c(y.id))
        
        // the number of distinct cj for j != i  
        val K = clusterVector.size
          
        val yCluster = getCluster(c(y.id))
        // innovation (new clusters)
        addAuxiliaryParameters(H)
        // assign y 
        c += (y.id -> max(calculateProbaWithLog(y)))

        if(c(y.id) >= K){
          val Beta = new BetaDistribution(1, gamma)
          val b = Beta.sample
          var cluster = getCluster(c(y.id))
          cluster.id = K
          cluster.size = 1
          cluster.beta = b * betaU
          c += (y.id -> K)
          addCluster(K, cluster)
          
          betaU = (1-b) * betaU
          
        }
        else
                   
          getCluster(c(y.id)).addPoint()
        
        for(j <- 0 to H - 1){
          
          val index = clusterVector.length - 1
          clusterVector.remove(index)
          
        }
        
      }

        alpha = alphaInference(alpha, 1, 0.5, clusterVector.size, nbPoint)
     
    }
     
    //calculate means of clusters
    for(cluster <- clusterVector)    { 
        
      var dataOfCluster = new ArrayBuffer[Point]
      // data of cluster
      for(d <- data)
        if (c(d.id) == cluster.id)
          dataOfCluster += d
        
      cluster.calculateMean(dataOfCluster)

    }
    
    var result = new ArrayBuffer[ResumeCluster]
    
    for(cluster <- clusterVector){
      
      if(cluster.size != 0){
      
      val sub = new ResumeCluster(workerId, cluster.id, cluster.size, cluster.y_ , cluster.phi)
      result += sub
      }
      
    }
    
    result

  }
  
  def addAuxiliaryParameters(H : Int) = {
    
    
    val N2 = new MultivariateNormalDistribution(m, sigma2)
    
    for (j <- 0 to H-1){

      val auxiliaryPhi = N2.sample()

      val id = clusterVector.length

      val a = new Array[Double](dim)
      for(i <- 0 to dim - 1)
        a(i) = 0

      val cluster = new Cluster(id, 0, a, auxiliaryPhi, betaU/H, true)

      clusterVector += cluster
    
    }
    
  }
  
  def calculateProbaWithLog(y : Point)={
    
    var proba : Map[Int, Double] = Map()
    for(cluster <- clusterVector){

        proba += (cluster.id -> (Math.log(cluster.size + alpha * cluster.beta) + Math.log(normalLikelihood(y.vector, cluster.phi, sigma1))- Math.log(nbPoint - 1 + alpha)))

    }
  
  val max = proba.values.max

    for ((cluster, p) <- proba){
      proba += cluster -> Math.exp(p-max)
  
    }
    
    val sum = proba.foldLeft(0.0){ case (a, (k, v)) => a + v } 
    
    for ((cluster, p) <- proba){
      proba += (cluster -> p/sum)
  
    }

    proba
  }
  
  def normalLikelihood(data : Array[Double], mean : Array[Double], sigma1 : Array[Array[Double]]) : Double ={
    var vraisemblance : Double = 1
    val N = new MultivariateNormalDistribution(mean, sigma1)
    N.density(data)
  
  }
  
  def max(proba : Map[Int, Double])={
  
    val max = proba.maxBy(_._2)._1
    
    max
  }
 
  def alphaInference(alpha : Double, a : Double, b : Double, k : Int, n : Int) = {
    val Beta = new BetaDistribution(alpha + 1, n)
    val eta = Beta.sample
    val pi = (a + k - 1)/(a + k - 1 + n * (b - Math.log(eta)))
   
    val rand = scala.util.Random
    val u = rand.nextDouble
   
   
    if (u <= pi){
      val Gamma = new GammaDistribution(a + k, b - Math.log(eta))
      Gamma.sample
    }else {
      val Gamma = new GammaDistribution(a + k - 1, b - Math.log(eta))
      Gamma.sample
    }
  }
  
  def contingencyTable(nbGlobalClusters : Int, nbRealClusters : Int) = {
     
    var mat = Array.fill[Int](nbGlobalClusters, nbRealClusters)(0)        

    for(y <- data)
      mat(c(y.id))(y.realCluster - 1) += 1
    
    mat  
        
  }

  def dataCluster(clusterId : Int) = {
    
    var dataOfCluster = new ArrayBuffer[Point]
    
    for(d <- data)
      if (c(d.id) == clusterId)
        dataOfCluster += d
      
    dataOfCluster
    
  }
 
 def localRSS() = {
   
    var RSS = 0.0
   
    for(clust <- clusterVector){
      
      val dataOfCluster = dataCluster(clust.id)
      
      for(y <- dataOfCluster)
        for(i <- 0 to dim - 1)
          RSS += Math.pow(y.vector(i) - clust.phi(i), 2) / E
      
    }
    
    RSS
    
 }
 
}
