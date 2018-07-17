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









There are a few ways to move along at a constant speed along a path whose "segments" are not a constant length - and it's not trivial to make them that way.

I would approach the problem by making a "Mover" of sorts which follows a "Path".

public interface Path<T> {
   public T getPoint(float delta);
}
public class Mover<T> {
   public Path<T> path;
   public T position;
   public T lookahead;
   public float lookaheadDelta;
   public float delta;
   public float speed;
   public void init() { ... }
   public void update(float dt) { ... }
}

(sorry if you haven't learned about generics/templating, just replace T with Vector3)

Most "path" implementations like Catmull-Rom spline can give you a point given a delta value that ranges from 0 to 1. Based off of this the Mover moves along the path looking "lookaheadDelta" ahead and moving towards that point. The smaller this value the smoother the movement, but if it's too small you might make to many path calculations per update.

The update method is important, it tries to move speed*dt units along the path. If the distance between lookahead and the current position is less than this value than you need to iteratively calculate the new lookahead based on the lookaheadDelta and continually move towards it.

It would look something like this:

public void update(float dt) {
  float move = dt * speed; // units to move
  while (move > 0.0) {
    // between current position and target
    float room = distance( lookahead, position ); 
    // how much we're actually moving this iteration
    float actual = Math.min( move, room );
    // the normalized vector between lookahead and position
    float direction = (lookahead - position) / room; 
    // move my position accordingly
    position += direction * actual;
    // update move to be the remaining amount we need to move
    move -= actual;
    // reset target
    if (actual != room) {
      delta += lookaheadDelta;
      lookahead = path.getPoint(delta);
    }
  }
}

And don't forget about initializing it!

public void init() {  
  position = path.getPoint(0.0);
  delta = lookaheadDelta;
  lookahead = path.getPoint(delta);
}

Don't forget to add in logic to check when delta >= 1.0 which means you're at the end!

You could even calculate lookaheadDelta based on the estimated length of the spline and the units of your game. I would guess 100 data points would be enough making 0.01 a good value for lookaheadDelta.

You could also cache "direction" and "room" after each lookahead calculation.

