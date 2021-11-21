#include <iostream>
#include <windows.h>
#include <vector>
#include <algorithm>
#include <queue>
#include <omp.h>

#define TRAIN_DATA 50000

using namespace std;

double start_time, end_time;

int* readData(const char* dataset_filename, int* N) {
    FILE* f;
    errno_t err = fopen_s(&f, "car_evaluation.csv", "r");
    if (err) {
        printf_s("The file car_evaluation was not opened\n");
        return 0;
    }
    fscanf_s(f, "%ld", N);
    int* data_points = (int*)malloc((*N) * 7 * sizeof(int));
    for (int i = 0; i < (*N) * 7; i++) {
        int val;
        fscanf_s(f, "%d", &val);
        *(data_points + i) = val;
    }

    fclose(f);
    return data_points;
}

double entropy(vector<int> freq, int sum) {
    double res = 0.0;
    for (int v : freq) {
        double prob = (v * 1.0) / (sum * 1.0);
        res -= prob * log(prob);
    }
    return res;
}

int getAt(int* data, int i, int j) {
    try
    {
        return *(data + i * 7 + j);
    }
    catch (const std::exception&)
    {
        printf("error att %d,%d\n", i, j);

    }
}

void uniqueValue(int* data, vector<int> ids, int attr, vector<int>& res) {
    res.clear();
    for (int id : ids) {
        res.push_back(getAt(data, id, attr));
    }
    sort(res.begin(), res.end());
    res.erase(unique(res.begin(), res.end()), res.end());
}

class TreeNode {
public:
    vector<int> ids;
    double entropy;
    int depth;
    vector<TreeNode> children;
    int splitAttribute;
    vector<int> order;
    int label;

    TreeNode() {}

    TreeNode(vector<int> ids, double entropy, int depth) {
        this->ids = ids;
        this->entropy = entropy;
        this->depth = depth;
        this->splitAttribute = -1;
        this->label = -1;
    }

    TreeNode(vector<int> ids, vector<TreeNode> children, double entropy, int depth) {
        this->ids = ids;           // index of data in this node
        this->entropy = entropy;   // entropy, will fill later
        this->depth = depth;       // distance to root node
        this->splitAttribute = -1; // which attribute is chosen, it non-leaf
        this->children = children; // list of its child nodes
        this->label = -1;       // label of node if it is a leaf
    }

    void setProperties(int splitAttribute, vector<int> order) {
        this->splitAttribute = splitAttribute;
        this->order = order;
    }

    void setLabel(int label) {
        this->label = label;
    }
};

class DecisionTreeID3 {
public:
    vector<int> ids;
    int maxDepth;
    int minSamplesSplit;
    double minGain;
    TreeNode root;
    int* data;
    vector<int> attributes;
    vector<int> labels;

    DecisionTreeID3(int maxDepth, int minSamplesSplit, double minGain) {
        this->maxDepth = maxDepth;
        this->minSamplesSplit = minSamplesSplit;
        this->minGain = minGain;
    }

    void fit(int* data_points) {
        this->data = data_points;
        this->attributes.insert(this->attributes.end(), { 0, 1, 2, 3, 4, 5 });
        this->labels.insert(this->labels.end(), { 1, 2, 3, 4 });
        for (int i = 0; i < TRAIN_DATA; i++) this->ids.push_back(i);
        this->root = TreeNode(ids, this->_entropy(ids), 0);
        queue<TreeNode> treeQueue; treeQueue.push(this->root);
        while (!treeQueue.empty()) {
            TreeNode node = treeQueue.front(); treeQueue.pop();
            //printf("%f \n", node.entropy);
            if (node.depth < this->maxDepth || node.entropy < this->minGain) {
                node.children = this->_split(node);
                if (node.children.size() == 0) { // leaf node
                    this->_set_label(node);
                }
                for (TreeNode treeNode : node.children)
                    treeQueue.push(treeNode);
            }
            else {
                this->_set_label(node);
            }
        }
    }

    double _entropy(vector<int> ids) {
        if (ids.size() == 0) return 0;
        vector<int> freq;
        int valueCount[6] = { 0, 0, 0, 0, 0, 0 };
        for (int id : ids) {
            valueCount[getAt(this->data, id, 6)]++;
        }
        int sum = 0;
        for (int i = 0; i < 6; i++) {
            if (valueCount[i] == 0) continue;
            freq.push_back(valueCount[i]);
            sum += valueCount[i];
        };
        return entropy(freq, sum);
    }

    void _set_label(TreeNode node) {
        int valueCount[6] = { 0, 0, 0, 0, 0, 0 };
        int maxPos = 0;
        for (int id : ids) {
            int val = getAt(this->data, id, 6);
            valueCount[val]++;
            if (valueCount[val] > valueCount[maxPos]) maxPos = val;
        }
        node.setLabel(maxPos); // most frequent label
    }

    vector<TreeNode> _split(TreeNode node) {
        vector<int> ids = node.ids;
        double bestGain = 0;
        int bestAttribute = -1;
        vector<int> order;
        vector<vector<int> > bestSplits;

        end_time = (double)clock() / (double)CLOCKS_PER_SEC;

        double executed_time = end_time - start_time;
        printf("Start find attr at %f \n", executed_time);

        omp_set_num_threads(6);
        #pragma omp parallel for
        for (int i = 0; i < this->attributes.size(); i++) {
            vector<int> subData;
            copy(ids.begin(), ids.end(), back_inserter(subData));
            int att = this->attributes.at(i);
            vector<int> values;
            uniqueValue(this->data, subData, i, values);
            if (values.size() == 1) continue; // entropy = 0
            vector<vector<int>> splits;
            splits.clear();
            for (int value : values) {
                vector<int> subIds;
                for (int j = 0; j < subData.size(); j++) {
                    if (getAt(this->data, subData.at(j), att) == value) {
                        subIds.push_back(subData.at(j));
                    }
                }
                splits.push_back(subIds);
            }

            // don't split if a node has too small number of points
            bool check = true;
            for (vector<int> split : splits) {
                if (split.size() < this->minSamplesSplit) {
                    check = false; break;
                }
            }
            if (!check) continue;

            // information gain
            double HxS = 0;
            for (vector<int> split : splits) {
                HxS += split.size() * this->_entropy(split) / ids.size();
            }
            double gain = node.entropy - HxS;
            if (gain < this->minGain) continue; // stop if small gain
            if (gain > bestGain) {
                bestGain = gain;
                bestSplits = splits;
                bestAttribute = att;
                order = values;
            }
            //Sleep(30);
        }
        end_time = (double)clock() / (double)CLOCKS_PER_SEC;

        executed_time = end_time - start_time;
        printf("Find best gain after %f return %f \n", executed_time, bestGain);
        node.setProperties(bestAttribute, order);
        vector<TreeNode> childNodes;
        for (vector<int> split : bestSplits) {
            childNodes.push_back(TreeNode(split, this->_entropy(split), node.depth + 1));
        }

        return childNodes;
    }
};

void decisionTree(int N, int* data_points) {
    DecisionTreeID3 tree = DecisionTreeID3(5, 2, 1e-4);
    tree.fit(data_points);
}

int main(int argc, char const* argv[]) {

    //---------------------------------------------------------------------//
    int N;
    const char* dataset_filename = "car_evaluation.csv";//argv[2];
    int* data_points;
    double executed_time;
    //---------------------------------------------------------------------//
    data_points = readData(dataset_filename, &N);

    printf("----------------Open MP Decision Tree Started------------------\n");
    start_time = (double)clock() / (double)CLOCKS_PER_SEC;
    decisionTree(N, data_points);
    end_time = (double)clock() / (double)CLOCKS_PER_SEC;

    executed_time = end_time - start_time;
    printf("Executed Time: %lf \n", executed_time);
    printf("-----------------------------Finished-----------------------------\n");
    return 0;
}