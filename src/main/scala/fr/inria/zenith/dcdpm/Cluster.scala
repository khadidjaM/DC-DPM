package fr.inria.zenith.dcdpm

import scala.collection.mutable.ArrayBuffer
import org.apache.commons.math3.distribution.MultivariateNormalDistribution
import Jama.Matrix

abstract class AbstractCluster(var id : Int, var size : Double, var y_ : Array[Double], var phi : Array[Double]) extends Serializable{
  
  var dim = y_.length
  
  def empty() = {size == 0}
  
  def calculateMean()={}
  
  def updateId(newId : Int) = {
    
      id = newId
      
  }
 
}

class Cluster(id : Int, s : Double, y : Array[Double], p : Array[Double], var beta : Double, var newCluster : Boolean) extends AbstractCluster(id, s, y , p){
 
  def isNew() = newCluster
    
  def addPoint() = {
      
    size += 1
      
  }
    
  def removePoint(pId : Int) = {
      
    size -= 1
      
  }
  
  def calculateMean(dataCluster : ArrayBuffer[Point])={
    
    if (empty)
      
      phi = Array.fill(dim)(0.0)
    
    else{
      
      val m = new Array[Double](dim)

      for(j <- 0 to dim - 1){

        var s : Double = 0

        for(d <- dataCluster)
            s = s + d.vector(j)

        m(j) = s / dataCluster.length

      }

      y_ = m
    }
    
  }
    
}

class GlobalCluster(var subIds : ArrayBuffer[(Int, Int)], id : Int, size : Double, y_ : Array[Double], p : Array[Double], var mean : Array[Double], var sigma2 : Array[Array[Double]])extends AbstractCluster(id, size, y_ , p){
  
  def updateY_(means : ArrayBuffer[_<:AbstractCluster])={
    
    var sum = Array.fill(dim)(0.0)

    var sumSize : Double = 0
    for (mean <- means){
      for(i <- 0 to dim - 1)
        sum(i) += mean.size * mean.y_(i)
      sumSize += mean.size
    }
    
    for(i <- 0 to dim - 1)
      y_(i) = sum(i) / sumSize
    
  }
  
  def updatePhi(sigma1 : Array[Array[Double]]) = {
   
  val sigma = postSigma(sigma1)
  
  postMean(sigma1, sigma)
   
  val N = new MultivariateNormalDistribution(mean, sigma)
  
    phi = N.sample
    
    sigma2 = sigma
      
  }
  
  def postMean(sigma1 : Array[Array[Double]], sigma : Array[Array[Double]]) = {
    
    val y = Array.ofDim[Double](dim, 1)
    for(i <- 0 to dim - 1)
      y(i)(0) = y_(i)
    

    var s1 = matrixProduct(inverseMatrix(sigma1), y)
    
    for (i <- 0 to dim - 1)
        s1(i)(0) *= size
      
    val m = Array.ofDim[Double](dim, 1)
    for(i <- 0 to dim - 1)
      m(i)(0) = mean(i)
    
    val s2 = matrixProduct(inverseMatrix(sigma2), m)
    
    val s = matrixProduct(sigma, matrixSum(s1, s2))    
    
    for(i <- 0 to dim - 1)
      mean(i) = s(i)(0)
  
  }
    
  def postSigma(sigma1 : Array[Array[Double]]) = {
    
    var sigma = inverseMatrix(sigma1)
    
    for (i <- 0 to dim - 1)
      for (j <- 0 to dim - 1)
        sigma(i)(j) *= size
    
    val s = inverseMatrix(matrixSum(sigma, inverseMatrix(sigma2)))
    
    s
  }
  
  def matrixProduct(A : Array[Array[Double]], B : Array[Array[Double]]) : Array[Array[Double]] = {
    
    val m = A.length
    val n = B.length
    val p = B(0).length
    var C = Array.ofDim[Double](m, p)
    
    for(i <- 0 to m -1)
      for(j <- 0 to p - 1)
        {
          var som : Double = 0
          
          for(k <- 0 to n - 1)
            som += A(i)(k) * B(k)(j)
          
          C(i)(j) = som
        }
    C
  }

  def matrixSum(A : Array[Array[Double]], B : Array[Array[Double]]) : Array[Array[Double]] = {
    
    val n = A.length
    val m = A(0).length
    var C = Array.ofDim[Double](n, m)

    for(i <- 0 to n - 1)
      for(j <- 0 to m - 1)
        C(i)(j) = A(i)(j) + B(i)(j)
        
    C 
  }
  
  def inverseMatrix (A : Array[Array[Double]]) = {
    
    val M = new Matrix(A)
    
    M.inverse.getArray
    
  }
  
  
}

class ResumeCluster(val workerId : Int, id : Int, size : Double, y_ : Array[Double], phi : Array[Double]) extends AbstractCluster(id, size, y_ , phi ){}
© 2018 GitHub,
