/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description:  main function to test the new read threading module implemented
 *    		    using template classes
 *
 *        Version:  1.0
 *        Created:  07/27/2015 12:58:38
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  JJ Chai (jjchai), jjchai01@gmail.com
 *   Organization:  ORNL
 *
 * =====================================================================================
 */


#include "Config.h"
#include "logcpp/log.h"
#include "tRead.h"
#include "namespace_jj.h"
#include "tDataSet.h"
#include "tEdge.h"
#include "tOverlapGraph.h"
#include "Utils.h"

int main() {
	CLOCKSTART;
	tDataSet<tRead<int> > * data = new tDataSet<tRead<int>>;
	for(int offset = 0; offset < 20; offset++){
		vector<int> v;
		for(int i = 0; i < 10; i++){
			v.push_back(((i-offset > 0 ) ? (i-offset):(offset-i)));
		}
		tRead<int> r(v);
		data->addRead(r);
	}
	for(int offset = 0; offset < 20; offset++){
		vector<int> v;
		for(int i = 0; i < 10; i++){
			v.push_back(((i-offset > 0 ) ? (i-offset):(offset-i)));
		}
		tRead<int> r(v);
		data->addRead(r);
	}

	data->rmDupReads();
	cout << "reads in the data set are sorted? " << std::boolalpha << data->isSorted() << endl;
	cout << "data does not have duplicate reads? " << std::boolalpha << data->isDupRemoved() << endl;
	cout << *data << endl;
	

	tOverlapGraph<int> * graph = new tOverlapGraph<int>(data, 1);
	vector<tEdge<tRead<int>>*> contigEdges;
	graph->getEdges(contigEdges);
	graph->printGraph("tmp.gdl", contigEdges);
	delete graph;
	delete data;

	CLOCKSTOP;
	return 0;
}
