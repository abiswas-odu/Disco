package tax;

import java.io.PrintStream;
import java.util.HashSet;
import java.util.List;
import java.util.regex.Pattern;

import stream.Read;

import fileIO.ReadWrite;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date Nov 30, 2015
 *
 */
public class TaxFilter {
	
	public static void main(String[] args){
		String regex=args[0];
		String s=args[1];
		Pattern regexPattern=(regex==null ? null : Pattern.compile(regex));
		boolean b=regexPattern.matcher(s).matches();
		System.err.println(regex);
		System.err.println(s);
		System.err.println(b);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public static TaxFilter makeFilter(String[] args){
		
		String names=null;
		String ids=null;

		String tableFile=null;
		String treeFile=null;

		int taxLevel=0;
		int reqLevel=0;
		boolean include=false;
		String regex=null;
		String contains=null;
		
		for(int i=0; i<args.length; i++){
			String arg=args[i];

			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);} //Strip leading hyphens
			
			if(a.equals("table") || a.equals("gi")){
				tableFile=b;
				if("auto".equalsIgnoreCase(b)){tableFile=TaxTree.defaultTableFile();}
			}else if(a.equals("tree") || a.equals("taxtree")){
				treeFile=b;
				if("auto".equalsIgnoreCase(b)){treeFile=TaxTree.defaultTreeFile();}
			}else if(a.equals("level") || a.equals("taxlevel")){
				if(Character.isDigit(b.charAt(0))){
					taxLevel=Integer.parseInt(b);
				}else{
					taxLevel=TaxTree.stringToLevel(b.toLowerCase());
				}
			}else if(a.equals("reqlevel") || a.equals("requiredlevel") || a.equals("reqlevels") || a.equals("requiredlevels")){
				String[] split2=b.toLowerCase().split(",");
				reqLevel=0;
				for(String s : split2){
					int level;
					if(Character.isDigit(s.charAt(0))){
						level=Integer.parseInt(s);
					}else{
						level=TaxTree.stringToLevel(s);
					}
					reqLevel|=(1<<level);
				}
			}else if(a.equals("name") || a.equals("names")){
				names=b;
			}else if(a.equals("regex")){
				regex=b;
			}else if(a.equals("contains")){
				contains=b;
			}else if(a.equals("include")){
				include=Tools.parseBoolean(b);
			}else if(a.equals("exclude")){
				include=!Tools.parseBoolean(b);
			}else if(a.equals("requirepresent")){
				REQUIRE_PRESENT=Tools.parseBoolean(b);
				TaxTree.SHOW_WARNINGS=REQUIRE_PRESENT;
			}else if(a.equals("id") || a.equals("ids") || a.equals("taxid") || a.equals("taxids")){
				ids=b;
			}
		}
		
		TaxFilter filter=new TaxFilter(tableFile, treeFile, taxLevel, reqLevel, include, null, regex, contains);
		filter.addNames(names);
		filter.addNumbers(ids);
		
		return filter;
	}
	
	/**
	 * Constructor.
	 * @param tree_
	 * @param taxLevel_
	 * @param include_
	 * @param taxSet_
	 */
	public TaxFilter(TaxTree tree_, int taxLevel_, int reqLevel_, boolean include_, HashSet<Integer> taxSet_, String regex_, String contains_){
		tree=tree_;
		taxLevel=taxLevel_;
		reqLevels=reqLevel_;
		include=include_;
		taxSet=(taxSet_==null ? new HashSet<Integer>() : taxSet_);
		regex=regex_;
		regexPattern=(regex==null ? null : Pattern.compile(regex));
		containsString=contains_;
	}
	
	/**
	 * Constructor.
	 * @param tableFile
	 * @param treeFile
	 * @param taxLevel_
	 * @param include_
	 * @param taxSet_
	 */
	public TaxFilter(String tableFile, String treeFile, int taxLevel_, int reqLevel_, boolean include_, HashSet<Integer> taxSet_, String regex_, String contains_){
		taxLevel=taxLevel_;
		reqLevels=reqLevel_;
		include=include_;
		taxSet=(taxSet_==null ? new HashSet<Integer>() : taxSet_);
		
		loadGiTable(tableFile);
		tree=loadTree(treeFile);
		regex=regex_;
		regexPattern=(regex==null ? null : Pattern.compile(regex));
		containsString=contains_;
	}
	
	public static boolean validArgument(String a){
		if(a.equals("table") || a.equals("gi")){
		}else if(a.equals("tree") || a.equals("taxtree")){
		}else if(a.equals("level") || a.equals("taxlevel")){
		}else if(a.equals("name") || a.equals("names")){
		}else if(a.equals("regex")){
		}else if(a.equals("contains")){
		}else if(a.equals("include")){
		}else if(a.equals("exclude")){
		}else if(a.equals("id") || a.equals("ids") || a.equals("taxid") || a.equals("taxids")){
		}else if(a.equals("requirepresent")){
		}else if(a.equals("reqlevel") || a.equals("requiredlevel") || a.equals("reqlevels") || a.equals("requiredlevels")){
		}else{
			return false;
		}
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	static void loadGiTable(String fname){
		if(fname==null){return;}
		if(PRINT_STUFF){outstream.println("Loading gi table.");}
		GiToNcbi.initialize(fname);
	}
	
	static TaxTree loadTree(String fname){
		if(fname==null){return null;}
		if(PRINT_STUFF){outstream.println("Loading tree.");}
		TaxTree tt=ReadWrite.read(TaxTree.class, fname, true);
		if(tt.nameMap==null){
			if(PRINT_STUFF){outstream.println("Hashing names.");}
			tt.hashNames();
		}
		assert(tt.nameMap!=null);
		return tt;
	}
	
	public void addNames(String names){
		if(names==null){return;}
		String[] array=names.split(",");
		for(String name : array){
			addName(name);
		}
	}
	
	public boolean addName(String name){
		{
			TaxNode tn=tree.getNode(name);
			if(tn!=null){return addNode(tn);}
		}
		List<TaxNode> list=tree.getNodesByName(name);
		boolean success=false;
		assert(list!=null) : "Could not find a node for '"+name+"'";
		for(TaxNode tn : list){
			success=addNode(tn)|success;
		}
		return success;
	}
	
	public void addNumbers(String numbers){
		if(numbers==null){return;}
		String[] array=numbers.split(",");
		for(String s : array){
			final int x=Integer.parseInt(s);
			addNumber(x);
		}
	}
	
	public boolean addNumber(int taxID){
		TaxNode tn=tree.getNode(taxID);
		assert(tn!=null) : "Could not find a node for '"+taxID+"'";
		return addNode(tn);
	}
	
	public boolean addNode(TaxNode tn){
		if(tn==null || tn.level>taxLevel){return false;}
		taxSet.add(tn.id);
		while(tn.id!=tn.pid && tn.level<=taxLevel){
			System.err.println("Added node "+tn);
			taxSet.add(tn.id);
			tn=tree.getNode(tn.pid);
		}
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Outer Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	public boolean passesFilter(final Read r){
		return passesFilter(r.id);
	}
	
	public boolean passesFilter(final String name){
		if(regexPattern!=null){
			boolean b=matchesRegex(name);
			if(b!=include){return false;}
		}
		if(containsString!=null){
			boolean b=containsString(name);
			if(b!=include){return false;}
		}
		if(taxSet.isEmpty() && reqLevels==0){return !include;}
		TaxNode tn=tree.getNode(name);
		if(tn==null){tn=tree.getNodeByName(name);}
		assert(tn!=null || !REQUIRE_PRESENT) : "Could not find node for '"+name+"'";
//		assert(false) : passesFilter(tn);
		return passesFilter(tn);
	}
	
	public boolean passesFilter(final int id){
		if(taxSet.isEmpty() && reqLevels==0){return !include;}
		TaxNode tn=tree.getNode(id);
		assert(tn!=null || !REQUIRE_PRESENT) : "Could not find node number "+id;
		return passesFilter(tn);
	}
	
	boolean passesFilter(final TaxNode tn0){
		TaxNode tn=tn0;
		if(taxSet.isEmpty() && reqLevels==0){return !include;}
		if(tn==null){
			assert(!REQUIRE_PRESENT) : "Null TaxNode.";
			return !include && reqLevels==0;
		}
		boolean found=taxSet.contains(tn.id);
//		System.err.println("found="+found+", node="+tn);
		int levels=1<<tn.level;
		while((!found || (levels&reqLevels)!=reqLevels) && tn.id!=tn.pid){
			tn=tree.getNode(tn.pid);
			levels|=(1<<tn.level);
			found=found||taxSet.contains(tn.id);
//			System.err.println("found="+found+", node="+tn);
		}
//		assert(false) : levels+", "+reqLevels+", "+tn0+", "+tree.getAncestors(tn0.pid);
		return include==found && (levels&reqLevels)==reqLevels;
	}
	
	boolean matchesRegex(String s){
		return regexPattern.matcher(s).matches();
	}
	
	boolean containsString(String s){
		return s.toLowerCase().contains(containsString);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public TaxTree tree(){return tree;}
	private final TaxTree tree;
	
	/** Level at which to filter */
	private final int taxLevel;
	
	/** Branch must contain ancestors at these levels (bitflag) */
	private final int reqLevels;
	
	/** Set of numeric NCBI TaxIDs */
	private final HashSet<Integer> taxSet;
	
	private final boolean include;
	
	private String regex;
	private final Pattern regexPattern;
	private final String containsString;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private static PrintStream outstream=System.err;
	
	/** Print loading messages */
	static boolean PRINT_STUFF=true;
	
	static boolean REQUIRE_PRESENT=true;

}
