**watershed algorithm for blob detection**

below we generate random 2d gaussians (left), detect peaks and use watershed for image segmentation (right)

<img src="https://user-images.githubusercontent.com/89974426/135768972-b4d92532-63a9-4311-aded-c31427ac2846.png" width=30% height=30%>    <img src="https://user-images.githubusercontent.com/89974426/135768998-61d06b70-86fe-4f41-a114-233e727871a5.png" width=30% height=30%>

**algorithm**
```
1.	Background noise estimation

•	Select smallest 1% pixels

•	Compute mean and variance from their 4-connectivity neighbors 

2.	Gaussian smoothing with estimated variance

3.	Convolution with Laplace filter, sharpening

4.	Replace pixels smaller than noise mean with 0

5.	Retrieve list of maximum regions with maximum filter from original picture

6.	Watershed from maximum regions
```


**addressing color invasion**

Illustrating color invasion

<img src="https://user-images.githubusercontent.com/89974426/135769771-4db3cea7-6ebe-49e9-b8d0-f96ab3642a14.png" width=30% height=30%> <img src="https://user-images.githubusercontent.com/89974426/135769776-845a1d11-0f5d-435b-bcf6-fb875a7e6d0b.png" width=30% height=30%> 

To address color invasion

•	Check before we label a region whether this region is in the 8-connectivity neighborhood of a region with this label. If not, we leave the region unlabeled. Other criteria for region labeling still hold: 
o	If a region has two higher neighbor region with different label, it is background
o	If a region has a higher neighbor which is background, it is background

•	Know when a blob is resolved. A blob is fully resolved when it is no more surrounded with unlabeled, in the 8-connectivity sense. As soon as a blob is fully resolved, we empty its queue. 

That's better

<img src="https://user-images.githubusercontent.com/89974426/135769831-1abab1ed-0673-4c2c-9cdb-35ba28357e7e.png" width=30% height=30%> 

But "coeluted zones" are detected poorly.

**hollow coeluted zones**

Another example 

<img src="https://user-images.githubusercontent.com/89974426/135770027-688f8c98-c889-4bae-80f5-6087eee1befd.png" width=70% height=70%>

coeluted zones are here

<img src="https://user-images.githubusercontent.com/89974426/135769968-576a431f-ebf2-4568-b0a7-930ee9a0515c.png" width=35% height=35%> 

or nicely seen in 3D


 <img src="https://user-images.githubusercontent.com/89974426/135770104-dd054f72-81c7-4d4f-b9f2-cc4df352d768.png" width=60% height=60%> 


Intuition

We have a contour for the coeluted zone. This is the ‘envelop for points already resolved as belonging to the blob’.

We wish to fill the zone with colors from the blobs in the zone. The demarcation line is the valley between adjacent blobs. 

We can distinguish between these two cases:

  •	Almost resolved zone, as zone with blobs 1 and 10 on figure 5b. Integration error for each of these blobs is approximately 15%. Valley could be defined as the set of lowest points from paths binding already resolved regions from one blob to already resolved regions from second blob.

  •	Poorly resolved zone as three blobs zone. 

   -Blobs 6 and 2 are said to be strongly coeluted. They could not be resolved separately, they are merged and the peaks apex with the highest volume only is kept

   -This conclusion holds at this scale. We don’t definitively discard the blob, but given preprocessing, this blob is not ‘seen’ by segmentation method.

Algorithm for coeluted zones filling
```
1.	Find coeluted zones

Blobs in a same coeluted zone are identified

2.	Find a hull

Quickhull algorithm, [7], computes convex hull for the coeluted zone. 
The convex hull is the smallest convex set containing all points already identified as belonging to blobs contained in this coeluted zone. 
Implementation details are in appendix.

Improvement could resort to alpha-shapes to closer fit the blob’s shape. 

3.	Find points in valley

4.	Complete valley 

From isolated points located in valley, complete valley finding the minimum spanning tree of these points

5.	Fill with color with valley as a separation
```

Convex hull

We try to join points in valleys in the most natural way. We build the minimum spanning tree of points in valley with Kruskal algorithm and we add to valleys all points on edges. 

<img src="https://user-images.githubusercontent.com/89974426/135770349-6f4cb9a2-196a-4749-add8-107e819997e9.png" width=30% height=30%>

Green dots fill in valleys smoothly.
