/*
Задача 5. Вариант 1. Минимальное остовное дерево

Дан неориентированный связный граф. Требуется найти вес минимального остовного дерева в этом графе.
Вариант 1. С помощью алгоритма Прима.
Вариант 2. С помощью алгоритма Крускала.
Вариант 3. С помощью алгоритма Борувки.
Ваш номер варианта прописан в ведомости.
*/

#include <iostream>
#include <vector>
#include <cassert>
#include <utility>
#include <set>
#include <tuple>
#include <climits>  // для INT_MAX
#include <random>
#include <algorithm>
#include <cmath>

struct IGraph {
	virtual ~IGraph() {}

	virtual void AddEdge(int from, int to, double weight) = 0;

	virtual int VerticesCount() const = 0;

	virtual std::vector<std::pair<double, int>> GetNextVerticesWithWeights(int vertex) const = 0;
	virtual std::vector<std::pair<double, int>> GetPrevVerticesWithWeights(int vertex) const = 0;
};

struct ListGraph : IGraph {
public:
	ListGraph(int size);
	ListGraph(const IGraph& graph);

	~ListGraph();

	void AddEdge(int from, int to, double weight) override;

	int VerticesCount() const override;

	void printGraph() const;

	std::vector<std::pair<double, int>> GetNextVerticesWithWeights(int vertex) const override;
	std::vector<std::pair<double, int>> GetPrevVerticesWithWeights(int vertex) const override;

private:
	std::vector<std::vector<std::pair<double, int>>> adjacencyLists;
};

ListGraph::ListGraph(int size) : adjacencyLists(size) {}

ListGraph::~ListGraph() {}

ListGraph::ListGraph(const IGraph& graph) : adjacencyLists(graph.VerticesCount()) {
	for (int i = 0; i < graph.VerticesCount(); ++i) {
		adjacencyLists[i] = graph.GetNextVerticesWithWeights(i);
	}
}

void ListGraph::AddEdge(int from, int to, double weight) {
	assert(from >= 0 && from < adjacencyLists.size());
	assert(to >= 0 && to < adjacencyLists.size());

	adjacencyLists[from].push_back(std::make_pair(weight, to));
	adjacencyLists[to].push_back(std::make_pair(weight, from));
}

int ListGraph::VerticesCount() const {
	return static_cast<int>(adjacencyLists.size());
}

std::vector<std::pair<double, int>> ListGraph::GetNextVerticesWithWeights(int vertex) const {
	assert(vertex >= 0 && vertex < adjacencyLists.size());
	return adjacencyLists[vertex];
}

std::vector<std::pair<double, int>> ListGraph::GetPrevVerticesWithWeights(int vertex) const {
	assert(vertex >= 0 && vertex < adjacencyLists.size());
	std::vector<std::pair<double, int>> prevVertices;
	for (int from = 0; from < adjacencyLists.size(); ++from) {
		for (int to = 0; to < adjacencyLists[from].size(); ++to) {
			if (adjacencyLists[from][to].second == vertex) {
				prevVertices.push_back(std::make_pair(adjacencyLists[from][to].first, from));
			}
		}
	}
	return prevVertices;
}

void ListGraph::printGraph() const {
	for (int v = 0; v < adjacencyLists.size(); ++v) {
		std::cout << "[" << v << "]: ";
		for (int u = 0; u < adjacencyLists[v].size(); ++u) {
			std::cout << "(" << adjacencyLists[v][u].second << ", " << adjacencyLists[v][u].first << ")";
		}
		std::cout << std::endl;
	}
}

bool relax(int u, int v, double weight_uv, std::vector<double>& min_e, std::vector<int>& p) {
	if (min_e[v] > weight_uv) {
		min_e[v] = weight_uv;
		p[v] = u;
		return true;
	}
	return false;
}

std::vector<int> MST_Prim(const ListGraph& graph, int s) {
	std::vector<double> min_e(graph.VerticesCount(), INT_MAX);
	std::vector<bool> visited(graph.VerticesCount(), false);
	std::vector<int> p(graph.VerticesCount(), -1);
	min_e[s] = 0;
	visited[s] = true;
	std::set<std::pair<double, int>> q;
	q.insert(std::make_pair(min_e[s], s));
	while (!q.empty()) {
		int u = (*(q.begin())).second;
		q.erase(*(q.begin()));
		for (auto nextVertex : graph.GetNextVerticesWithWeights(u)) {
			int v = nextVertex.second;
			if (visited[v]) {
				continue;
			}
			double w_uv = nextVertex.first;
			double min_e_v_old = min_e[v];
			if (min_e[v] == INT_MAX) {
				min_e[v] = w_uv;
				p[v] = u;
				q.insert(std::make_pair(min_e[v], v));
			}
			else if (relax(u, v, w_uv, min_e, p)) {
				q.erase(std::make_pair(min_e_v_old, v));
				q.insert(std::make_pair(min_e[v], v));
			}
		}
		visited[u] = true;
	}

	return p;
}

double computeRouteLength(const std::vector<std::pair<double, double>>& vertices, const std::vector<int>& p) {
	double length = 0;
	for (int i = 0; i < p.size() - 1; ++i) {
		length += sqrt(pow(vertices[p[i]].first - vertices[p[i + 1]].first, 2) + pow(vertices[p[i]].second - vertices[p[i + 1]].second, 2));
	}
	length += sqrt(pow(vertices[p[p.size() - 1]].first - vertices[p[0]].first, 2) + pow(vertices[p[p.size() - 1]].second - vertices[p[0]].second, 2)); // Замыкание маршрута
	return length;
}

double bruteForce(const std::vector<std::pair<double, double>>& vertices) {
	std::vector<int> permutation;
	for (int i = 0; i < vertices.size(); ++i) {
		permutation.push_back(i);
	}
	double min_dist = computeRouteLength(vertices, permutation);
	do {
		double cur_dist = computeRouteLength(vertices, permutation);
		if (cur_dist < min_dist) {
			min_dist = cur_dist;
		}
	} while (std::next_permutation(permutation.begin(), permutation.end()));

	return min_dist;
}

void dfs(int node, const std::vector<std::vector<int>>& tree, std::vector<int>& preorder) {
	preorder.push_back(node);

	for (int child : tree[node]) {
		dfs(child, tree, preorder);
	}
}

double approximateMST(const std::vector<std::pair<double, double>>& vertices) {
	int start = 0;
	int n = vertices.size();
	ListGraph listGraph(n);

	for (int i = 0; i < n; ++i) {
		for (int j = i + 1; j < n; ++j) {
			listGraph.AddEdge(i, j, sqrt(pow(vertices[i].first - vertices[j].first, 2) + pow(vertices[i].second - vertices[j].second, 2)));
		}
	}

	std::vector<int> p = MST_Prim(listGraph, start);
	std::vector<std::vector<int>> tree(n);

	for (int i = 0; i < n; ++i) {
		if (p[i] != -1) {
			tree[p[i]].push_back(i);
		}
	}

	std::vector<int> preorder;

	int root = -1;
	for (int i = 0; i < n; ++i) {
		if (p[i] == -1) {
			root = i;
			break;
		}
	}

	if (root != -1) {
		dfs(root, tree, preorder);
	}

	return computeRouteLength(vertices, preorder);
}

int main()
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dist(-1.0, 1.0);

	int n_iters = 25;

	for (int n = 4; n < 10; ++n) {
		double m1 = 0;
		double m2 = 0;
		double s = 0;

		std::cout << "-------------------" << std::endl;
		std::cout << "n = " << n << std::endl << std::endl;

		for (int iter = 0; iter < n_iters; ++iter) {
			std::vector<std::pair<double, double>> vertices;
			for (int i = 0; i < n; ++i) {
				double u, v;
				double s = 0.0;
				while (s == 0 || s > 1) {
					u = dist(gen);
					v = dist(gen);
					s = u * u + v * v;
				}
				double z0 = u * sqrt(-2 * log(s) / s);
				double z1 = v * sqrt(-2 * log(s) / s);

				vertices.push_back(std::make_pair(z0, z1));
			}

			double C_star = bruteForce(vertices);
			double C = approximateMST(vertices);
			std::cout << "Iter " << iter << ":" << std::endl;
			std::cout << "Brute force: " << C_star << std::endl;
			std::cout << "Approximate: " << C << std::endl;
			std::cout << std::endl;

			double app_coef = C / C_star > C_star / C ? C / C_star : C_star / C;

			m1 += app_coef;
			m2 += app_coef * app_coef;
		}

		m1 /= n_iters;
		m2 /= n_iters;

		s = sqrt(m2 - m1 * m1); // M(X^2) - M(X)^2 - дисперсия

		std::cout << "Mean value of approximation coefficient for n = " << n << ": " << m1 << std::endl;
		std::cout << "Standart deviation of approximation coefficient for n = " << n << ": " << s << std::endl;
		std::cout << std::endl;
	}
}