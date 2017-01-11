package stream;

import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;

public class ReadStreamStringWriter extends ReadStreamWriter {

	@Deprecated
	public ReadStreamStringWriter(String fname_, boolean read1_, int bufferSize, boolean allowSubprocess_){
		this(FileFormat.testOutput(fname_, FileFormat.BREAD, null, allowSubprocess_, false, false, true), null, read1_, bufferSize, null, false);
	}
	
	public ReadStreamStringWriter(FileFormat ff, String qfname_, boolean read1_, int bufferSize, CharSequence header, boolean useSharedHeader){
		super(ff, qfname_, read1_, bufferSize, header, true, true, useSharedHeader);
	}

	@Override
	public void run() {
		
		if(!OUTPUT_SAM && !OUTPUT_FASTQ && !OUTPUT_FASTA && !OUTPUT_ATTACHMENT && !OUTPUT_HEADER && !OUTPUT_ONELINE){
			if(OUTPUT_INTERLEAVED){
//				assert(false) : OUTPUT_SAM+", "+OUTPUT_FASTQ+", "+OUTPUT_FASTA+", "+OUTPUT_ATTACHMENT+", "+OUTPUT_INTERLEAVED+", "+SITES_ONLY;
				myWriter.print("#INTERLEAVED\n");
			}
			if(SITES_ONLY){
				myWriter.println("#"+SiteScore.header());
			}else if(!OUTPUT_ATTACHMENT){
				myWriter.println("#"+Read.header());
			}
		}
		
		Job job=null;
		while(job==null){
			try {
				job=queue.take();
//				job.list=queue.take();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		while(job!=null && !job.poison){
//			System.err.println("Processing job "+job);
			if(!job.isEmpty()){
				
				if(myQWriter!=null){
					if(read1){
						for(final Read r : job.list){
							if(r!=null){
								{
									CharSequence cs=(r.bases==null ? "\n" : toQualitySB(r.quality, r.length(), FASTA_WRAP).append('\n'));
									myQWriter.print('>');
									myQWriter.println(r.id);
									myQWriter.print(cs);
								}
								Read r2=r.mate;
								if(OUTPUT_INTERLEAVED && r2!=null){
									CharSequence cs=(r2.bases==null ? "\n" : toQualitySB(r2.quality, r2.length(), FASTA_WRAP).append('\n'));
									myQWriter.print('>');
									myQWriter.println(r2.id);
									myQWriter.print(cs);
								}
							}
						}
					}else{
						for(final Read r1 : job.list){
							if(r1!=null){
								final Read r2=r1.mate;
								assert(r2!=null && r2.mate==r1 && r2!=r1) : r1.toText(false);
								CharSequence cs=(r2.bases==null ? "\n" : toQualitySB(r2.quality, r2.length(), FASTA_WRAP).append('\n'));
								myQWriter.print('>');
								myQWriter.println(r2.id);
								myQWriter.print(cs);
							}
						}
					}
				}
//				assert(false) : OUTPUT_SAM+", "+SITES_ONLY+", "+OUTPUT_FASTQ+", "+OUTPUT_FASTA+", "+OUTPUT_ATTACHMENT+"\n"+job.list.get(0).obj+"\n"+job.list.get(0);
				if(OUTPUT_SAM){
					assert(read1);
					for(final Read r : job.list){
						Read r2=(r==null ? null : r.mate);

						SamLine sl1=(r==null ? null : (USE_ATTACHED_SAMLINE && r.obj!=null ? (SamLine)r.obj : new SamLine(r, 0)));
						SamLine sl2=(r2==null ? null : (USE_ATTACHED_SAMLINE && r2.obj!=null ? (SamLine)r2.obj : new SamLine(r2, 1)));
						
						if(r!=null){
							job.writer.print(sl1.toText().append('\n'));

							readsWritten++;
							basesWritten+=(r.bases!=null ? r.length() : 0);
							ArrayList<SiteScore> list=r.sites;
							if(OUTPUT_SAM_SECONDARY_ALIGNMENTS && list!=null && list.size()>1){
								final Read clone=r.clone();
								for(int i=1; i<list.size(); i++){
									SiteScore ss=list.get(i);
									clone.match=null;
									clone.setFromSite(ss);
									clone.setSecondary(true);
//									assert(false) : r.mapScore+"\n"+ss.header()+"\n"+r.list+"\n";
									SamLine sl=new SamLine(clone, 0);
									assert(!sl.primary());
//									sl.setPrimary(false);
									
									job.writer.print(sl.toText().append('\n'));

//									readsWritten++;
//									basesWritten+=(r.bases!=null ? r.length() : 0);
//									validReadsWritten+=(r.valid() && r.mapped() ? 1 : 0);
//									validBasesWritten+=(r.valid() && r.mapped() && r.bases!=null ? r.length() : 0);
								}
							}
						}
						if(r2!=null){
							if(!SamLine.KEEP_NAMES && sl1!=null && ((sl2.qname==null) || !sl2.qname.equals(sl1.qname))){
								sl2.qname=sl1.qname;
							}
							job.writer.print(sl2.toText().append('\n'));

							readsWritten++;
							basesWritten+=(r2.bases!=null ? r2.length() : 0);
							
							ArrayList<SiteScore> list=r2.sites;
							if(OUTPUT_SAM_SECONDARY_ALIGNMENTS && list!=null && list.size()>1){
								final Read clone=r2.clone();
								for(int i=1; i<list.size(); i++){
									SiteScore ss=list.get(i);
									clone.match=null;
									clone.setFromSite(ss);
									clone.setSecondary(true);
//									assert(false) : r.mapScore+"\n"+ss.header()+"\n"+r.list+"\n";
									SamLine sl=new SamLine(clone, 0);
									assert(!sl.primary());
//									sl.setPrimary(false);
									if(!SamLine.KEEP_NAMES && sl1!=null && ((sl2.qname==null) || !sl2.qname.equals(sl1.qname))){
										sl2.qname=sl1.qname;
									}
									job.writer.print(sl.toText().append('\n'));

//									readsWritten++;
//									basesWritten+=(r.bases!=null ? r.length() : 0);
//									validReadsWritten+=(r.valid() && r.mapped() ? 1 : 0);
//									validBasesWritten+=(r.valid() && r.mapped() && r.bases!=null ? r.length() : 0);
								}
							}
						}
						
					}
				}else if(SITES_ONLY){
					assert(read1);
					for(final Read r : job.list){
						Read r2=(r==null ? null : r.mate);
						
						if(r!=null && r.sites!=null){
							job.writer.print(r.toSites().append('\n'));

							readsWritten++;
							basesWritten+=(r.bases!=null ? r.length() : 0);
						}
						if(r2!=null){
							job.writer.print(r2.toSites().append('\n'));

							readsWritten++;
							basesWritten+=(r2.bases!=null ? r2.length() : 0);
						}
					}
				}else if(OUTPUT_FASTQ){
					if(read1){
						for(final Read r : job.list){
							if(r!=null){
								job.writer.print(r.toFastq().append('\n'));
								readsWritten++;
								basesWritten+=(r.bases!=null ? r.length() : 0);
								Read r2=r.mate;
								if(OUTPUT_INTERLEAVED && r2!=null){
									job.writer.print(r2.toFastq().append('\n'));
									readsWritten++;
									basesWritten+=(r2.bases!=null ? r2.length() : 0);
								}
							}
						}
					}else{
						for(final Read r1 : job.list){
							if(r1!=null){
								final Read r2=r1.mate;
								assert(r2!=null && r2.mate==r1 && r2!=r1) : r1.toText(false);
								job.writer.print(r2.toFastq().append('\n'));
								readsWritten++;
								basesWritten+=(r2.bases!=null ? r2.length() : 0);
							}
						}
					}
				}else if(OUTPUT_FASTA){
					if(read1){
						for(final Read r : job.list){
							if(r!=null){
								job.writer.print(r.toFasta(FASTA_WRAP).append('\n'));
								readsWritten++;
								basesWritten+=(r.bases!=null ? r.length() : 0);
								Read r2=r.mate;
								if(OUTPUT_INTERLEAVED && r2!=null){
									job.writer.print(r2.toFasta(FASTA_WRAP).append('\n'));
									readsWritten++;
									basesWritten+=(r2.bases!=null ? r2.length() : 0);
								}
							}
						}
					}else{
						for(final Read r1 : job.list){
							if(r1!=null){
								final Read r2=r1.mate;
								assert(r2!=null && r2.mate==r1 && r2!=r1) : r1.toText(false);
								job.writer.print(r2.toFasta(FASTA_WRAP).append('\n'));
								readsWritten++;
								basesWritten+=(r2.bases!=null ? r2.length() : 0);
							}
						}
					}
				}else if(OUTPUT_ONELINE){
					if(read1){
						for(final Read r : job.list){
							if(r!=null){
								job.writer.print(r.id);
								job.writer.print('\t');
								job.writer.print(new String(r.bases));
								job.writer.print('\n');
								readsWritten++;
								basesWritten+=(r.bases!=null ? r.length() : 0);
								Read r2=r.mate;
								if(OUTPUT_INTERLEAVED && r2!=null){
									job.writer.print(r2.id);
									job.writer.print('\t');
									job.writer.print(new String(r2.bases));
									job.writer.print('\n');
									readsWritten++;
									basesWritten+=(r2.bases!=null ? r2.length() : 0);
								}
							}
						}
					}else{
						for(final Read r1 : job.list){
							if(r1!=null){
								final Read r2=r1.mate;
								assert(r2!=null && r2.mate==r1 && r2!=r1) : r1.toText(false);
								job.writer.print(r2.id);
								job.writer.print('\t');
								job.writer.print(new String(r2.bases));
								job.writer.print('\n');
								readsWritten++;
								basesWritten+=(r2.bases!=null ? r2.length() : 0);
							}
						}
					}
				}else if(OUTPUT_ATTACHMENT){
					if(read1){
						for(final Read r : job.list){
							if(r!=null){
								job.writer.println(r.obj==null ? "." : r.obj.toString());
								readsWritten++;
								Read r2=r.mate;
								if(OUTPUT_INTERLEAVED && r2!=null){
									job.writer.println(r2.obj==null ? "." : r2.obj.toString());
									readsWritten++;
								}
							}
						}
					}else{
						for(final Read r1 : job.list){
							if(r1!=null){
								final Read r2=r1.mate;
//								assert(r2!=null && r2.mate==r1 && r2!=r1) : r1.toText(false);
								if(r2!=null){
									job.writer.println(r2.obj==null ? "." : r2.obj.toString());
									readsWritten++;
								}else{
									job.writer.println(".");
								}
							}
						}
					}
				}else if(OUTPUT_HEADER){
					if(read1){
						for(final Read r : job.list){
							if(r!=null){
								job.writer.println(r.id);
								readsWritten++;
								Read r2=r.mate;
								if(OUTPUT_INTERLEAVED && r2!=null){
									job.writer.println(r.id);
									readsWritten++;
								}
							}
						}
					}else{
						for(final Read r1 : job.list){
							if(r1!=null){
								final Read r2=r1.mate;
//								assert(r2!=null && r2.mate==r1 && r2!=r1) : r1.toText(false);
								if(r2!=null){
									job.writer.println(r2.id);
									readsWritten++;
								}else{
									job.writer.println(".");
								}
							}
						}
					}
				}else{
					if(read1){
						for(final Read r : job.list){
							if(r!=null){
								job.writer.print(r.toText(true).append('\n'));
								readsWritten++;
								basesWritten+=(r.bases!=null ? r.length() : 0);
								Read r2=r.mate;
								if(OUTPUT_INTERLEAVED && r2!=null){
									job.writer.print(r2.toText(true).append('\n'));
									readsWritten++;
									basesWritten+=(r2.bases!=null ? r2.length() : 0);
								}
								
							}
						}
					}else{
						for(final Read r1 : job.list){
							if(r1!=null){
								final Read r2=r1.mate;
//								assert(r2!=null && r2.mate==r1 && r2!=r1) : r1.toText(false);
								if(r2!=null){
									job.writer.print(r2.toText(true).append('\n'));
									readsWritten++;
									basesWritten+=(r2.bases!=null ? r2.length() : 0);
								}else{
									job.writer.print(".\n");
								}
							}
						}
					}
				}
			}
			if(job.close){
				assert(job.writer!=null && job.writer!=myWriter);
				ReadWrite.finishWriting(job.writer, job.outstream, fname, allowSubprocess);
			}
			
			job=null;
			while(job==null){
				try {
					job=queue.take();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		
		if(myWriter!=null){
			ReadWrite.finishWriting(myWriter, myOutstream, fname, allowSubprocess);
		}
		if(myQWriter!=null){
			ReadWrite.finishWriting(myQWriter, myQOutstream, qfname, allowSubprocess);
		}
		finishedSuccessfully=true;
	}
	
}
