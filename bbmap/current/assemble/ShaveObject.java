package assemble;

import java.io.PrintStream;

/**
 * Holds constants for shaving.
 * @author Brian Bushnell
 * @date Jul 20, 2015
 *
 */
public abstract class ShaveObject {
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print messages to this stream */
	static PrintStream outstream=System.err;
	
	/** Explore codes */
	public static final int KEEP_GOING=0, DEAD_END=1, TOO_SHORT=2, TOO_LONG=3, TOO_DEEP=4, FORWARD_BRANCH=5, BACKWARD_BRANCH=6, LOOP=7;
	public static final int STATUS_UNEXPLORED=0, STATUS_EXPLORED=1, STATUS_REMOVE=2, STATUS_KEEP=3;
	
	public static boolean printEventCounts=false;
	
	/** Verbose messages */
	public static boolean verbose=false;
	/** Debugging verbose messages */
	public static boolean verbose2=false;
	
}
