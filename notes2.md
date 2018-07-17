Brief description of the code:
1) The 3rd parameter is wantedDistance which is the distance between two nodes.
2) You can see the const DIVISOR - it is just to increase the accuracy of the distance because the algo I tried is this:-
a) fineGrainDistance = wantedDistance/DIVISOR
b) Calculate the actual point (Vector3) on the Spline where the distance = fineGrainDistance
c) Do this for DIVISOR times.. you get the actual node which has distance 'wantedDistance'
d) Repeat until reach the end of spline.

So if your spline is curvy (lots of sharp bends) and the desired distance between node is big, set DIVISOR to big value.

void getEqualDistanceSpline(SimpleSpline& splineSrc, SimpleSpline* splineDest, Real wantedDistance)
{
	Real lastInterpPoint = 0.0;
	Real length;
	Vector3 start = splineSrc.getPoint(0);
	Vector3 end;
	Real wantedDistanceSquared = wantedDistance*wantedDistance;

	splineDest->setAutoCalculate(false);
	splineDest->addPoint(start);

	for (int j = 1; j < splineSrc.getNumPoints();) {

		// first find the points where the length exceed wanted length..
		end = splineSrc.getPoint(j);
		length = (end-start).squaredLength();

		while (length < wantedDistanceSquared && j < splineSrc.getNumPoints()-1) {
			end = splineSrc.getPoint(++j);
			length = (end-start).squaredLength();
			// if enter the loops then we have to reset lastInterPoint..
			lastInterpPoint = 0.0;
		}

      // Moved to the end of the "for" loop and changed (see below)
//		if (j == splineSrc.getNumPoints() -1)
//			break;

		// okay found it.. lets refine
		Real partStart = lastInterpPoint;
		Real partEnd = 1.0;
		Real partMid;
		Vector3 partPoint;
		Real partLen;
		const Vector3& refPoint = splineSrc.interpolate(j-1, lastInterpPoint);
		Real squaredDist = wantedDistance-(start-refPoint).length();
		squaredDist *= squaredDist;
		
		do {
			partMid = (partStart+partEnd)/2;
			partPoint = splineSrc.interpolate(j-1, partMid);
			partLen = (partPoint-refPoint).squaredLength();
			if (fabs(partLen-squaredDist)< 1 || fabs(partStart-partEnd) < 1e-5)
				break;
			if (partLen > squaredDist)
				partEnd = partMid;
			else
				partStart = partMid;
		} while (true);

		// once we reach here.. the exact point has been discovered..
		start = splineSrc.interpolate(j-1, partMid);
		//LOG("\tstart = " + StringConverter::toString(start) + ", lastInterpPoint = " + StringConverter::toString(partMid));
		// and remember the last interpolation point
		lastInterpPoint = partMid;

		splineDest->addPoint(start);

      //
      // Moved from above; this was exiting too soon and 
      // not including nodes for the last leg of a spline.
      //
      if ((length < rDistBetweenNodesSquared) && (j == this->getNumPoints() - 1))
	   {
         break;
	   }
	}
	splineDest->recalcTangents();
}
