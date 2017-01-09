package assemble;

public class ErrorTracker {
	
	public ErrorTracker(){
		
	}
	
	public void clear(){
		clearDetected();
		clearCorrected();
		
		rollback=false;
		suspected=0;
		marked=0;
	}
	
	public void clearDetected(){
		detectedPincer=0;
		detectedTail=0;
		detectedBrute=0;
		detectedReassemble=0;
	}

	public void clearCorrected() {
		correctedPincer=0;
		correctedTail=0;
		correctedBrute=0;
		correctedReassembleInner=0;
		correctedReassembleOuter=0;
	}
	
	public int corrected(){
		return correctedPincer+correctedTail+correctedBrute+correctedReassembleInner+correctedReassembleOuter;
	}
	
	public int detected(){
		return detectedPincer+correctedTail+detectedReassemble;
	}
	
	public int correctedReassemble(){
		return correctedReassembleInner+correctedReassembleOuter;
	}
	
	public String toString(){
		StringBuilder sb=new StringBuilder();
		sb.append("suspected         \t").append(suspected).append('\n');
		sb.append("detectedPincer    \t").append(detectedPincer).append('\n');
		sb.append("detectedTail      \t").append(detectedTail).append('\n');
//		sb.append("detectedBrute     \t").append(detectedBrute).append('\n');
		sb.append("detectedReassemble\t").append(detectedReassemble).append('\n');
		sb.append("correctedPincer   \t").append(correctedPincer).append('\n');
		sb.append("correctedTail     \t").append(correctedTail).append('\n');
//		sb.append("correctedBrute    \t").append(correctedBrute).append('\n');
		sb.append("correctedReassembleInner\t").append(correctedReassembleInner).append('\n');
		sb.append("correctedReassembleOuter\t").append(correctedReassembleOuter).append('\n');
		sb.append("marked            \t").append(marked);
		return sb.toString();
	}

	public int suspected;
	
	public int detectedPincer;
	public int detectedTail;
	public int detectedBrute;
	public int detectedReassemble;
	
	public int correctedPincer;
	public int correctedTail;
	public int correctedBrute;
	public int correctedReassembleInner;
	public int correctedReassembleOuter;
	
	public int marked;
	
	public boolean rollback=false;
	
}
