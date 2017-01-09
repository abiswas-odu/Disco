package hiseq;

import shared.Tools;

public class FlowcellCoordinate {
	
	public FlowcellCoordinate() {}
	
	public FlowcellCoordinate(String id) {
		setFrom(id);
	}
	
	public float distance(FlowcellCoordinate fc){
		assert(isSet());
		assert(fc.isSet());
		
		if(lane!=fc.lane){return big;}
		
		long a=Tools.absdif(x, fc.x), b=Tools.absdif(y, fc.y);
		if(tile!=fc.tile){
			return spanTiles ? Tools.min(a, b) : big;
		}
		return (float)Math.sqrt(a*a+b*b);
		
		//Hard to say...  could consider adjacent tiles?
		//TODO: Ensure coordinates are not continuous across tiles.
//		if(tile!=fc.tile){
//			if(allowAdjacentTiles && Tools.absdif(tile, fc.tile)<2){return Tools.min(x-fc.x, y-fc.y);}
//			return big;
//		}
//		
//		long a=x-fc.x, b=y-fc.y;
//		return (float)Math.sqrt(a*a+b*b);
	}

	//2402:6:1101:6337:2237/1
	//MISEQ08:172:000000000-ABYD0:1:1101:18147:1925 1:N:0:TGGATATGCGCCAATT
	//HISEQ07:419:HBFNEADXX:1:1101:1238:2072
	public void setFrom(String id){
		final int lim=id.length();
		
		int i=0;
		int current=0;
		while(i<lim && id.charAt(i)!=' ' && id.charAt(i)!='/'){i++;}
		for(int semis=0; i>=0; i--){
			if(id.charAt(i)==':'){
				semis++;
				if(semis==4){break;}
			}
		}
		i++;
		
		assert(Character.isDigit(id.charAt(i))) : id;
		while(i<lim && Character.isDigit(id.charAt(i))){
			current=current*10+(id.charAt(i)-'0');
			i++;
		}
		lane=current;
		current=0;
		i++;
		
		if(!Character.isDigit(id.charAt(i))){//Hiseq 3000?
			while(i<lim && id.charAt(i)!=':'){i++;}
			i++;

			assert(Character.isDigit(id.charAt(i))) : id;
			while(i<lim && Character.isDigit(id.charAt(i))){
				current=current*10+(id.charAt(i)-'0');
				i++;
			}
			lane=current;
			current=0;
			i++;
		}

		assert(Character.isDigit(id.charAt(i))) : id;
		while(i<lim && Character.isDigit(id.charAt(i))){
			current=current*10+(id.charAt(i)-'0');
			i++;
		}
		tile=current;
		current=0;
		i++;

		assert(Character.isDigit(id.charAt(i))) : id;
		while(i<lim && Character.isDigit(id.charAt(i))){
			current=current*10+(id.charAt(i)-'0');
			i++;
		}
		x=current;
		current=0;
		i++;

		assert(Character.isDigit(id.charAt(i))) : id;
		while(i<lim && Character.isDigit(id.charAt(i))){
			current=current*10+(id.charAt(i)-'0');
			i++;
		}
		y=current;
		current=0;
		i++;
	}
	
	public boolean isSet(){
		return lane>=0 && tile>=0 && x>=0 && y>=0;
	}
	
	public int lane=-1;
	public int tile=-1;
	public int x=-1;
	public int y=-1;
	
	public static final float big=10000000;
	public static boolean spanTiles=true;
	
}
