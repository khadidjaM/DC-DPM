/*
 * Copyright 2018 Khadidja Meguelati <khadidja.meguelati@inria.fr>.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package fr.inria.zenith.dcdpm

import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._
import org.apache.spark.SparkConf 
import scala.collection.mutable.ArrayBuffer
import org.apache.commons.math3.distribution.GammaDistribution
import org.apache.spark.mllib.clustering.KMeans
import org.apache.spark.mllib.linalg.Vectors

object App { 
  
  def convertIteratorToArrayBuffer(data : Iterator[Point]) = {
    
    var dataBuffer = new ArrayBuffer[Point]
    var index = 0
    
    while (data.hasNext){
      
      dataBuffer.insert(index, data.next)
      index += 1
      
    }
      
    dataBuffer
    
  }
  
  def arrayToArrayBuffer(A : ArrayBuffer[ResumeCluster])= {
    
    var B = new ArrayBuffer[ResumeCluster]
    
    for (i <- 0 to A.length - 1)
      if(A(i).size > 0){
        B += A(i)
      
      }
    
    B
     
  }
 
  def toPoint(d : Array[Double], i : Int, r : Int)={
    
    val p = new Point(i, d, r)
    
    p
  }
  
  def dirichletSample(alpha : Array[Double]) = {
    
    var x = new Array[Double](alpha.length)
    var y = new Array[Double](alpha.length)
    
    for(i <- 0 to alpha.length - 1){
      
      val gamma = new GammaDistribution(alpha(i), 1)
      y(i) = gamma.sample
      
    }
    
    val sum = y.foldLeft(0.0){ case (a, b) => a + b } 
    
    for(i <- 0 to alpha.length - 1)
      x(i) = y(i)/sum
    
    x
      
  }
 
  def lineData(line : String, dim : Int, withReal : Int) = {
    
    val d = new Array[Double](dim)
    val temp = line.split(" ")
    
    for (j <- 0 to dim - 1)
      d(j) = temp(j).toDouble
    if(withReal == 1 )
      (d, temp(dim).toInt)
    else
      (d, 0)
    
  }
 
  def unionBuffers( A: Array[ArrayBuffer[ResumeCluster]]) = {
      
    var B = new ArrayBuffer[ResumeCluster]
      
    for(a <- A)
      B ++= a
        
    B
      
   }
  
  def main(args: Array[String]):Unit={
    
    if(args.length < 2) {
            System.err.println( "Use: DCDPM <variance in clusters> <variance between centers> <dimensions> <number of workers> <number of distributions> <target to data file> <number of clusters for Kmeans> <number of real clusters> <real clusters are known>")
            System.exit(1)
    }
      
    val conf = new SparkConf().setAppName("ParallelDP")
    val sc = new SparkContext(conf) 

    val E = args(0).toDouble
    val S = args(1).toDouble
    val dim = args(2).toInt
    val nbWorkers = args(3).toInt
    val nbDist = args(4).toInt
    val data = sc.textFile(args(5), nbWorkers).persist
    var numClusters = args(6).toInt
    val numRealClusters = args(7).toInt
    val withReal = args(8).toInt
    //val file = args(9)
    
    val n = data.count.toInt
    
    val parsedData = data.map(s => Vectors.dense(s.split(' ').take(dim).map(_.toDouble)))

    val numKmeansIterations = 10
    
    var gamma : Double = 1
    // intialize by Kmeans
    val model = KMeans.train(parsedData, numClusters, numKmeansIterations)
    
    val centers = model.clusterCenters
    
    val clusterInd = model.predict(parsedData)

    val clusterSizes = clusterInd.countByValue.values
     
    parsedData.unpersist()
    
    val dataRDD = data.zipWithIndex().map{ case (line, i) => {
          val p = lineData(line, dim, withReal) 
          toPoint(p._1, i.toInt, p._2)}}.persist
    
    val m = new Array[Double](dim)
    for(i <- 0 to dim - 1)
      m(i) = 0
  
    var sigma1 = Array.ofDim[Double](dim, dim)
    for(i <- 0 to dim - 1)
      for(j <- 0 to dim - 1)
        if(i == j)
          sigma1(i)(j)= E
        else
          sigma1(i)(j) = 0

    var sigma2 = Array.ofDim[Double](dim, dim)
    for(i <- 0 to dim - 1)
      for(j <- 0 to dim - 1)
        if(i == j)
          sigma2(i)(j)= S
        else
          sigma2(i)(j) = 0
    // create workers
    var workers = dataRDD.mapPartitionsWithIndex{(index, dataPartition) => {Iterator(new Worker(index, convertIteratorToArrayBuffer(dataPartition), gamma, E, S, 1))}}.persist
    
   // val writer = new PrintWriter(new File(file))
  //  val writer2 = new PrintWriter(new File(file+"C"))
    
    var time : Long = 0
    
    var beginTime = System.currentTimeMillis()
 
    var effectives = new Array[Double](numClusters + 1) 
    var  i = 0
    for (n <- clusterSizes){
      effectives(i) = n.toDouble
      i += 1
    }
    effectives(i) = gamma
    // sample Beta
    var Beta = dirichletSample(effectives)
    var betaU : Double = Beta(numClusters)

    var initial = new ArrayBuffer[(Array[Double], Double)]
    
      for(i <- 0 to centers.length-1){

        initial += centers(i).toArray -> Beta(i)
        
      }
    // initialze mean and sigma2 for each global cluster
    var parameters = new ArrayBuffer[(Array[Double], Array[Array[Double]])]
    for(i <-0 to numClusters - 1)
      parameters += m -> sigma2
      
    var initialClusters = new ArrayBuffer[GlobalCluster]
    
    var globalClusters = new ArrayBuffer[((Array[Double], Double), ArrayBuffer[(Int, Int)])]
    // initialize workers
    workers = workers.mapPartitions{ w => {Iterator(w.next.initializeWithClustering(initial, betaU))}}.persist
    
    for(j <- 0 to nbDist - 1){
      
      //val Alpha = workers.mapPartitions{ w => {Iterator(w.next.getAlpha)}}.collect
      
      if(j != 0){
        
        beginTime = System.currentTimeMillis()
       
        effectives = new Array[Double](numClusters + 1) 
        var i = 0
        for (c <- initialClusters){
          effectives(i) = c.size.toDouble
          i += 1
        }
        effectives(i) = gamma
        // sample Beta
        Beta = dirichletSample(effectives)
        betaU = Beta(numClusters)

        globalClusters = new ArrayBuffer[((Array[Double], Double), ArrayBuffer[(Int, Int)])]

        parameters = new ArrayBuffer[(Array[Double], Array[Array[Double]])]

        for(cluster <- initialClusters){

          globalClusters += cluster.phi -> Beta(cluster.id) -> cluster.subIds
          // post mean and sigma2 for each global cluster
          parameters += cluster.mean -> cluster.sigma2

        }
        // initialize workers with the new global clusters
        workers = workers.mapPartitions{ w => {Iterator(w.next.startWithClustering(globalClusters, betaU))}}.persist
        
        time += System.currentTimeMillis() - beginTime
     /*   
        writer.write("phis : ")
        for(clust <- initialClusters){           
          writer.write("\n")
          for(j <- 0 to dim - 1)
            writer.write(clust.phi(j)+ " ")
          writer.write("\n")
          
          writer.write("bar{y} : ")          
          for(j <- 0 to dim - 1)
            writer.write(clust.y_(j)+ " ")
          writer.write("\n")
      
      }
      
        */
        
        if(withReal == 1){
        
          val contingencyTables = workers.mapPartitions{ w => {Iterator(w.next.contingencyTable(globalClusters.length, numRealClusters))}}.collect

          val ARI = adjustedRandIndex(contingencyTables, nbWorkers, globalClusters.length, numRealClusters, n)

          println("ARI : " + ARI)

         // writer.write("ARI : " + ARI +" ")
        
        }
        
        val RSS = workers.mapPartitions{ w => {Iterator(w.next.localRSS)}}.collect.sum
        
        println("RSS : " + RSS)
        
       // writer.write("RSS : " + RSS)
        
      //  writer.write("\n")
        
      /*
        if(j == nbDist - 1){
          
          val c = workers.mapPartitions{ w => {Iterator(w.next.getC)}}.collect
          for(cj <- c)
            for(ci <- cj){
              writer2.write(ci._1 + " " + ci._2)
              writer2.write("\n")
            }
        
        }
        */
        
        beginTime = System.currentTimeMillis()
        
      }
      
      val resultClustering = workers.mapPartitions{ w => {Iterator(w.next.gibbsSampling())}}.collect
      
      val result = arrayToArrayBuffer(unionBuffers(resultClustering))
      
      var crp = new Master(E, S, n, gamma, parameters, dim)
      crp.updateExistingClasses(numClusters, result)
      crp.initialize      
      
       initialClusters.clear
       
       initialClusters = crp.gibbsSampling()
       
      numClusters = initialClusters.length
      
      gamma = crp.getGamma
       
      time += System.currentTimeMillis() - beginTime
             
      println("distribution : " + j + " nb Clusters : " + numClusters + " time : "  + time)
      /*
      writer.write("dist : " + j)
      writer.write("\n")
      writer.write(" alphas : ")
        for(a <- Alpha)
          writer.write(a + " ")
        writer.write("\n")
      
      writer.write("gamma : " + gamma)
      writer.write("\n")
      writer.write("nb Clusters : " + numClusters)
      writer.write("\n")
      writer.write("time : "  + time +" ")
      
      writer.write("\n")
      
      */
      
      time = 0
      
    }
    
  //  writer.close
  //  writer2.close
    
  }
 
  def Cn2(n : Int) : Double = n*(n-1)/2
  
  def adjustedRandIndex(contingencyTables : Array[Array[Array[Int]]], nbPart : Int, nbGlobalClusters : Int, numRealClusters : Int, n: Int)= {
    
    var globalContingencyTable = Array.fill[Int](nbGlobalClusters, numRealClusters)(0)
        
        for(m <- 0 to nbPart - 1)
          for(i <- 0 to nbGlobalClusters - 1)
            for(j <- 0 to numRealClusters - 1)
              globalContingencyTable(i)(j) += contingencyTables(m)(i)(j)
        
        var a  = Array.fill[Int](nbGlobalClusters)(0)
        var b  = Array.fill[Int](numRealClusters)(0)
        
        for(i <- 0 to nbGlobalClusters - 1)
          for(j <- 0 to numRealClusters - 1){
            a(i) += globalContingencyTable(i)(j)
            b(j) += globalContingencyTable(i)(j)
          }
        
        var index = 0.0
        var sumCai2 = 0.0
        var sumCbj2 = 0.0
        
        for(i <- 0 to nbGlobalClusters - 1){
          sumCai2 += Cn2(a(i))
          for(j <- 0 to numRealClusters - 1){
            index += Cn2(globalContingencyTable(i)(j))
            if(i == 0)
              sumCbj2 += Cn2(b(j))
          }
          
        }
        
    val expectedIndex = sumCai2 * sumCbj2 / Cn2(n)
        
    (index - expectedIndex) / (((sumCai2 + sumCbj2) / 2) - expectedIndex)
    
  }

}
