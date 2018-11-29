# DC-DPM

This is a Distributed Clustering with Spark based on Dirichlet Process Mixture, this approach is described in the following paper:

Khadidja Meguelati, Benedicte Fontez, Nadine Hilgert, and Florent Masseglia.2019. Dirichlet Process Mixture Models made Scalable and Effective by meansof Massive Distribution. InProceedings of ACM SAC Conference (SAC’19).

Please kindly cite our paper if the code helps you. Thank you.

## Requirements
DC-DPM works with [Apache Spark](http://spark.apache.org). In order to run it you must download and install [Spark Release 2.0.0](https://spark.apache.org/releases/spark-release-2-0-0.html).
The code is written in [Scala](https://www.scala-lang.org/), install [Scala 2.11.6](https://www.scala-lang.org/download/2.11.6.html)

## Building
We use maven to build it, Use the given [pom.xml](https://github.com/anonymeDoc/DC-DPM/blob/master/pom.xml) file to build an executable jar containing all the dependencies.

## Use
To execute DC-DPM use the following command :
```
$SPARK_HOME/bin/spark-submit --class "com.mycompany.dcdpm.App" DCDPM-jar-with-dependencies.jar <variance in clusters> <variance between centers> <dimensions> <number of workers> <number of distributions> <target to data file> <number of clusters for Kmeans> <number of real clusters> <real clusters are known>
```
### Necessary parameters
1. **variance in clusters:** We suppose that data are generated from a normal distribution, we need a covariance matrix with n dimensions which is an identity matrix with the value σ² in the diagonal. You should give the value of σ² 
2. **variance between centers:** We suppose that centers are generated from a normal distribution, we need a covariance matrix with n dimensions which is an identity matrix with the value σ² in the diagonal. You should give the value of σ² 
3. **dimensions:** the number of dimensions
5. **number of workers:**
6. **number of distributions:** in each distribution we perform several iterations of Gibbs Sampling on each worker and a synchronsation at the master level  
7. **target to data file:** The data file should be as follow :
  * each data in a line
  * values are seperated by space " "
  * if the ground truth is known, the data file should contain the label of the real cluster for each data in the last column, see [data with known ground truth.txt](https://github.com/anonymeDoc/DC-DPM/blob/master/data%20with%20known%20ground%20truth.txt) an example with 3 real clusters, and [data with unknown ground truth.txt](https://github.com/anonymeDoc/DC-DPM/blob/master/data%20with%20unknown%20ground%20truth.txt) an example of the other case.
8. **number of clusters for Kmeans:** the initialization Of DPM is done by a K-means step, you should indicate the number of clusters for Kmeans initialization
9. **number of real clusters:** if the ground truth is known, indicate the number of real clusters, else you can enter 0 
10. **real clusters are known:** if the ground truth is known, enter 1 else you can enter 0
