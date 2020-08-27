#include <iostream>
#include <vector>
#include <cassert>
#include <ctime>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <climits>
#include "PriorityQueue.h"
#include "Data.h"
#include "Graph.h"
#include "ShortestPath.h"

using namespace std;


const int INFINITY = INT_MAX;

vector<int> ShortestPath::path(const Graph &graph, int u, int w)
{
	vector<int> route;

	if (dijkstra(graph, u, w) == INFINITY)
	{
		return route;
	}
	
	int cur = w;
	do
	{
		route.push_back(cur);
		cur = traces[cur];

	} while (cur != -1);

	reverse(route.begin(), route.end());

	return route;
}

int ShortestPath::pathSize(const Graph &graph, int u, int w)
{
	return dijkstra(graph, u, w);
}

int ShortestPath::dijkstra(const Graph &graph, int u, int w)
{
	int size = graph.getVerticeNum();
	assert(u>=0 && u<size && w>=0 && w<size);

	// initiate traces;
	traces.clear();
	traces.resize(size, -1);

	PriorityQueue pQueue(graph.getVertices());
	pQueue.changePriority(QueueNode(u, 0));

	while (pQueue.size() > 0)
	{
		QueueNode top = pQueue.top();
		pQueue.pop();

		if (top.priority == INFINITY)
		{
			return INFINITY;
		}

		if (top.symbol == w)
		{
			return top.priority;
		}

		// relax path
		auto neighbors = graph.neighbors(top.symbol);
		QueueNode node;
		for (int i=0; i<neighbors.size(); i++)
		{
			node.symbol = neighbors[i];
			if (pQueue.contain(node))
			{
				int edge = graph.getEdgeValue(top.symbol, node.symbol);
				if (edge + top.priority < node.priority)
				{
					traces[node.symbol] = top.symbol;
					node.priority = edge + top.priority;
					pQueue.changePriority(node);
				}
			}
		}
	}

	return INFINITY;
}

PriorityQueue::PriorityQueue(const vector<int> &symbols)
{
	for (int i=0; i<symbols.size(); i++) 
	{
		if (indices.count(symbols[i]) == 0)
		{
			minHeap.push_back(QueueNode(symbols[i]));
			indices[symbols[i]] = i;
		}
	}

	minHeapfy();
}

void PriorityQueue::changePriority(const QueueNode &node)
{
	// if the node does not exist, return false;
	if (indices.count(node.symbol) == 0)
	{
		return;
	}

	int index = indices[node.symbol];
	if (minHeap[index].priority < node.priority)
	{
		minHeap[index].priority = node.priority;
		downwards(index);
	}
	else if (minHeap[index].priority > node.priority)
	{
		minHeap[index].priority = node.priority;
		upwards(index);
	}
}

void PriorityQueue::pop()
{
	if (size() <= 1)
	{
		minHeap.clear();
		indices.clear();
		return;
	}

	QueueNode &last = minHeap.back();
	indices.erase(minHeap[0].symbol);
	indices[last.symbol] = 0;
	minHeap[0].symbol = last.symbol;
	minHeap[0].priority = last.priority;
	minHeap.pop_back();

	downwards(0);
}

bool PriorityQueue::insert(const QueueNode &node)
{
	// if the node exists, return false;
	if (indices.count(node.symbol) > 0)
	{
		return false;
	}

	minHeap.push_back(node);

	int index = size() - 1;
	indices[node.symbol] = index;
	upwards(index);

	return true;
}

void PriorityQueue::minHeapfy()
{
	// nodes which have children need downwards heapfy
	int start = (size() - 2) / 2;

	for (int i=start; i>=0; i--)
	{
		downwards(i);
	}
}

void PriorityQueue::upwards(int i)
{
	if (i <= 0)
		return;

	int parent = (i - 1) / 2;

	if (minHeap[parent].priority > minHeap[i].priority)
	{
		swap(parent, i);
		upwards(parent);
	}
}

void PriorityQueue::downwards(int i)
{
	int left = i * 2 + 1;
	int right = i * 2 + 2;

	if (right < size())
	{
		if (minHeap[i].priority <= minHeap[left].priority &&
				minHeap[i].priority <= minHeap[right].priority)
		{
			return;
		}

		int smaller = 
			minHeap[left].priority<minHeap[right].priority ? left : right;

		swap(i, smaller);
		downwards(smaller);
	}
	else if (left < size())
	{
		if (minHeap[i].priority <= minHeap[left].priority)
		{
			return;
		}

		swap(i, left);
		downwards(left);
	}
}


Graph::Graph(int verticeNum): verticeNum(verticeNum)
{
	if (verticeNum <= 0)
	{
		this->verticeNum = 0;
		return;
	}

	adMatrix = vector<vector<int> >(verticeNum, vector<int>(verticeNum));
}

Graph::Graph(int verticeNum, double density):
	edgeNum(0), verticeNum(verticeNum)
{
	if (verticeNum <= 0)
	{
		this->verticeNum = 0;
		return;
	}

	adMatrix = vector<vector<int> >(verticeNum, vector<int>(verticeNum));

	srand(time(0));
	const int RANGE = 10;
	for (int i=0; i<verticeNum-1; i++)
	{
		for (int j=1; j<verticeNum; j++)
		{
			double prob = (rand() % 100) / 100.0;
			
			if (prob < density)
			{
				int value = rand() % RANGE + 1;
				addEdge(i, j, value);
			}
		}
	}
}

bool Graph::adjacent(int x, int y) const
{
	assert(x>=0 && x<verticeNum && y>=0 && y<verticeNum);

	return adMatrix[x][y] > 0;
}

vector<int> Graph::neighbors(int x) const
{
	assert(x>=0 && x<verticeNum);

	vector<int> list;
	for (int i=0; i<verticeNum; i++)
	{
		if (adMatrix[x][i] > 0)
		{
			list.push_back(i);
		}
	}

	return list;
}

vector<int> Graph::getVertices() const
{
	vector<int> vertices;

	for (int i=0; i<verticeNum; i++)
	{
		vertices.push_back(i);
	}

	return vertices;
}

bool Graph::addEdge(int x, int y, int value)
{
	assert(x>=0 && x<verticeNum && y>=0 && y<verticeNum);

	if (adMatrix[x][y] > 0)
		return false;

	adMatrix[x][y] = value;
	adMatrix[y][x] = value;

	return true;
}

bool Graph::deleteEdge(int x, int y)
{
	assert(x>=0 && x<verticeNum && y>=0 && y<verticeNum);

	if (adMatrix[x][y] > 0) 
	{
		adMatrix[x][y] = 0;
		adMatrix[y][x] = 0;
		return true;
	}

	return false;
}

int Graph::getEdgeValue(int x, int y) const
{
	assert(x>=0 && x<verticeNum && y>=0 && y<verticeNum);

	return adMatrix[x][y];
}

void Graph::setEdgeValue(int x, int y, int value)
{
	assert(x>=0 && x<verticeNum && y>=0 && y<verticeNum);

	adMatrix[x][y] = value;
	adMatrix[y][x] = value;
}


int main()
{
	int sum = 0;
	int count = 0;
	double densities[2] = {0.2, 0.4};
	int verticeNum = 50;
	ShortestPath sp;
	int i;
	int n;
	int pathSize;

	for (i=0; i<2; i++)
	{
		Graph graph(verticeNum, densities[i]);

		for (n=1; n<50; n++)
		{
			pathSize = sp.pathSize(graph, 0, n);
			if (pathSize != INFINITY)
			{
				sum += pathSize;
				count++;
			}
		}

		cout << "For the graph with density = " << densities[i];
		cout << ", the average path length = " << 
			(static_cast<double>(sum) / count) << endl;
	}

	return 0;
}
