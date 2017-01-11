package var2;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;

import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Tools;
import stream.Read;
import stream.ReadInputStream;
import stream.SamLine;

public class ScafMap {
	
	/*--------------------------------------------------------------*/
	/*----------------        Construction          ----------------*/
	/*--------------------------------------------------------------*/
	
	public ScafMap(){}

	public static ScafMap loadHeader(String fname){return loadHeader(fname, null);}
	public static ScafMap loadHeader(FileFormat ff){return loadHeader(ff, null);}
	
	public static ScafMap loadHeader(String fname, ScafMap scafMap){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.SAM, null, false, false);
		return loadHeader(ff, scafMap);
	}
	
	public static ScafMap loadHeader(FileFormat ff, ScafMap scafMap){
		ByteFile bf=ByteFile.makeByteFile(ff, false);
		if(scafMap==null){scafMap=new ScafMap();}
		byte[] line=bf.nextLine();
		while(line!=null && line.length>0){
			//assert(false) : new String(line);
			if(Tools.startsWith(line, "@SQ\t")){
				scafMap.add(line);
			}else if(line[0]!='@'){
				break;
			}
			line=bf.nextLine();
		}
		bf.close();

		return scafMap;
	}

	public static ScafMap loadReference(String fname){return loadReference(fname, null);}
	public static ScafMap loadReference(FileFormat ff){return loadReference(ff, null);}
	
	public static ScafMap loadReference(String fname, ScafMap scafMap){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
		return loadReference(ff, scafMap);
	}
	
	public static ScafMap loadReference(FileFormat ff, ScafMap map){
		ArrayList<Read> reads=ReadInputStream.toReads(ff, -1);
		if(map==null){map=new ScafMap();}
		for(Read r : reads){
			map.addScaffold(r);
		}
		return map;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Inherited           ----------------*/
	/*--------------------------------------------------------------*/
	
	public void clear(){
		map.clear();
		alt.clear();
		list.clear();
	}
	
	public int size(){return list.size();}
	
//	public boolean containsKey(String s){
//		return getScaffold(s)!=null;
//	}
	
	public Set<String> keySet() {
		return map.keySet();
	}
	
	public Set<String> altKeySet() {
		return alt.keySet();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Adders             ----------------*/
	/*--------------------------------------------------------------*/
	
	public Scaffold addScaffold(Read r){
		Scaffold scaf=map.get(r.id);
		if(scaf==null){
			scaf=new Scaffold(r.id, size(), r.length());
			add(scaf);
		}
		scaf.bases=r.bases;
		assert(scaf.bases.length==scaf.length);
		return scaf;
	}
	
	public Scaffold add(byte[] line){
		Scaffold scaf=new Scaffold(line, size());
		Scaffold old=map.get(scaf.name);
		if(old!=null){return old;}
		return add(scaf);
	}
	
	public Scaffold add(String s, int len){
		Scaffold scaf=map.get(s);
		if(scaf!=null){return scaf;}
		scaf=new Scaffold(s, size(), len);
		return add(scaf);
	}
	
	private Scaffold add(Scaffold scaf){
		assert(!map.containsKey(scaf.name));
		assert(size()==scaf.number);
		
		list.add(scaf);
		map.put(scaf.name, scaf);
		String s=scaf.name;
		if(TRIM_WHITESPACE_ALSO){
			for(int i=0; i<s.length(); i++){
				if(Character.isWhitespace(s.charAt(i))){
					String s2=s.substring(0, i);
					boolean b=alt.containsKey(s2);
					assert(!b);
					if(!b){
						alt.put(s2, scaf);
					}
					break;
				}
			}
		}
		return scaf;
	}	
	
	public void addCoverage(SamLine sl){
		if(!sl.mapped()){return;}
		Scaffold scaf=getScaffold(sl);
		scaf.add(sl);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/
	
	public int getNumber(String s){
		Scaffold value=getScaffold(s);
		return value==null ? -1 : value.number;
	}
	
	public String getName(int number){
		return number>=size() ? null : list.get(number).name;
	}
	
	public int getLength(int number){
		return number>=size() ? null : list.get(number).length;
	}
	
	public Scaffold getScaffold(int number){
		return number>=size() ? null : list.get(number);
	}
	
	public Scaffold getScaffold(String s){
		Scaffold value=map.get(s);
		if(value==null){value=alt.get(s);}
		if(value==null && TRIM_WHITESPACE_ALSO){
			int index=-1;
			for(int i=0; i<s.length(); i++){
				if(Character.isWhitespace(s.charAt(i))){
					index=i;
					break;
				}
			}
			if(index>0){
				String sub=s.substring(0, index);
				value=alt.get(sub);
				if(value==null){value=map.get(sub);}
			}
			assert(value!=null) : s+"\n"+index+"\n"+s.substring(0, Tools.max(1, index))
				+"\n"+keySet()+"\n"+altKeySet()+"\n";
		}
		assert(value!=null) : s+"\n"+keySet()+"\n"+altKeySet()+"\n";
		return value;
	}
	
	public Scaffold getScaffold(SamLine sl){
		return getScaffold(sl.rnameS());
	}
	
	public int getCoverage(Var v){
		Scaffold scaf=getScaffold(v.scafnum);
		return scaf.calcCoverage(v);
	}
	
	public long lengthSum() {
		long sum=0;
		for(Scaffold scaf : list){
			sum+=scaf.length;
		}
		return sum;
	}
	
	public void clearCoverage() {
		for(Scaffold scaf : list){scaf.clearCoverage();}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/

	final ArrayList<Scaffold> list=new ArrayList<Scaffold>();
	final HashMap<String, Scaffold> map=new HashMap<String, Scaffold>();
	private final HashMap<String, Scaffold> alt=new HashMap<String, Scaffold>();
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	public static ScafMap defaultScafMap=null;
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	public static boolean TRIM_WHITESPACE_ALSO=true;
	
}
