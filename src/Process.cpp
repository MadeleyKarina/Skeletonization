#include "Process.h"

Process::Process(){
    cout << "------------------------------------"<< endl;
    cout << "No initialization of current shape and skeleton type"<< endl;
}
Process::Process(string loadFile){
    currentShape = loadFile;
    cout << "constructor: " << currentShape <<endl;
}
Process::Process(char typeMatrix,  bool numVertDB){
    typeAdjMatrix = typeMatrix;
    changeNumberVert = numVertDB;
    currentShape = "/Dissim.off";
    cout << "constructor: " << typeAdjMatrix <<endl;
}
Process::Process(string loadFile, bool skelProperties, char typeMatrix){
    currentShape = loadFile;
    typeSkeleton = skelProperties;
    typeAdjMatrix = typeMatrix;
    changeNumberVert = false;
    currDiag = "corr_skel2.cgal";
    cout << "constructor: " << typeAdjMatrix <<endl;
}
Process::Process(string loadFile, bool skelProperties, char typeMatrix, bool numVertDB){
    currentShape = loadFile;
    typeSkeleton = skelProperties;
    typeAdjMatrix = typeMatrix;
    changeNumberVert = numVertDB;
    cout << "constructor: " << typeAdjMatrix <<endl;
}

//Euclidean distance
//template <typename T>
double	Process::vectors_distance(const std::vector<float>& a, const std::vector<float>& b)
{
    std::vector<float>	auxiliary;

    std::transform (a.begin(), a.end(), b.begin(), std::back_inserter(auxiliary),//
    [](float element1, float element2) {return pow((element1-element2),2);});
    auxiliary.shrink_to_fit();

    return  std::sqrt(std::accumulate(auxiliary.begin(), auxiliary.end(), 0.0));
} // end template vectors_distance
string Process::getShapeName(string loadFile){
    int start = loadFile.find_last_of('/');
    string shape= loadFile.substr(start+1);
    shape.erase(shape.find('.'),4);
    return shape;
}
string Process::getResultFileName(string loadFile){
    string shape= getShapeName(loadFile);
    cout << "Shape: "<< shape << endl;
    string nameFileMatrix;
    switch(typeAdjMatrix){
        case 'L': nameFileMatrix = "MatrixLaplaceLast/";
        break;
        case 'C': nameFileMatrix = "MatrixConnectivity/";
        break;
        case 'A': nameFileMatrix = "MatrixWeight/";
        break;
        case 'W': nameFileMatrix = "MatrixWeight_W/";
        break;
        default :
             cout << "Invalid type of matrix name" << endl;
    }
    //cout << nameFileMatrix << endl;
    std::string resultFile;
    if(changeNumberVert == true){
        resultFile = PATH_RESULT + string("Matrix_400_vertex/")
                + nameFileMatrix +shape+string("_")+to_string(numberEigen);
    }
    else{
        resultFile = PATH_RESULT + nameFileMatrix +shape+string("_")+to_string(numberEigen);
    }
    //cout <<resultFile<<endl;
    return resultFile;
}
const char * Process::getEigenFileName(string loadFile){
    std::string resultFile = getResultFileName(loadFile);
    string actualEigen = resultFile +"_eigen.cgal";
    cout << actualEigen << endl;
    const char *nameEigen = actualEigen.c_str();
    return nameEigen;
}
void Process::push_backMat(Eigen::MatrixXd& m, Vector3d&& values, std::size_t row)
{
    if(row >= m.rows()) {
        m.conservativeResize(row + 1, Eigen::NoChange);
    }
    m.row(row) = values;
}
vector<int> Process::getValency(Eigen::MatrixXd graph){
  vector<int> valency(graph.rows(),0);
  for (int i = 0; i < graph.rows(); i++){
      int n = 0;
      for (int j = 0; j < graph.rows(); j++){
          if(graph(i,j) == -1 )
          {
              n = n +1;
          }
      }
      valency[i]=n;
      //cout << n <<endl;
      //graph(i,i) = n;
  }
  return valency;
}
Eigen::MatrixXd Process::getLapGraph(Eigen::MatrixXd vertices, vector<int> edgesItem){
    int val;
    switch(typeAdjMatrix){
        case 'L': val = -1;
        break;
        case 'C': val = 1;
        break;
        case 'N': val = -1;
        break;
        case 'A': val = -1;
        break;
        case 'W': val = -1;
        break;
        default :
             cout << "Invalid type of matrix" << endl;
    }
    cout << "Graphhhhhhh "<<val << endl;
    Eigen::MatrixXd graph(vertices.rows() ,vertices.rows());
    graph= MatrixXd::Zero(vertices.rows() ,vertices.rows());
    for (int i = 0; i < edgesItem.size(); i+=2){
        graph(edgesItem[i],edgesItem[i+1])=val;
        graph(edgesItem[i+1],edgesItem[i])=val;
    }
    if(typeAdjMatrix=='A' || typeAdjMatrix == 'W'){
      vector<int> valency = getValency(graph);
      cout << "valency " << valency.size() <<endl;
      for (int i = 0; i < edgesItem.size(); i+=2){
          float weight = sqrt(valency[edgesItem[i]]*valency[edgesItem[i+1]]);
          //cout << "weight: "<< valency[edgesItem[i]]<<" "<<valency[edgesItem[i+1]]<<" "<<weight << endl;
          graph(edgesItem[i],edgesItem[i+1])=val/weight;
          graph(edgesItem[i+1],edgesItem[i])=val/weight;

          graph(i,i) = 1;

      }
      for(int i = 0; i < vertices.rows(); i++){
        if(typeAdjMatrix == 'W'){
          float aa = 1.0/valency[i];
          graph(i,i) = 1.0 - (1.0/valency[i]);
          //cout << 1.0 - (1.0/valency[i]) << " "<<i << "= i"<<endl;
        }
        else{
          graph(i,i) = 1;
        }
      }
    }
    if(typeAdjMatrix=='L'){
      vector<int> valency = getValency(graph);
      for(int i = 0; i < valency.size(); i++){
        graph(i,i) = valency[i];
      }
      // for (int i = 0; i < graph.rows(); i++){
      //       int n = 0;
      //       for (int j = 0; j < graph.rows(); j++){
      //           if(graph(i,j) == -1 )
      //           {
      //               n = n +1;
      //           }
      //       }
      //
      //       graph(i,i) = n;
      //   }
    }
    if(typeAdjMatrix=='N'){
        cout << "count Surf" <<endl;
//        vector<int> diag = countSurfVert(vertices,currDiag);
        vector<float> diag = distSurfVert(vertices,currDiag);
        vector<float>::iterator it = diag.begin();
        int i = 0;
        for (; it != diag.end(); it++ ){
           graph(i,i) = *it;
           i++;
        }
    }
    //std::cout << "Here is the matrix m:\n" << graph << std::endl;
    return graph;
}
Eigen::MatrixXd Process::readEdges(string fileCgaL){
    int cont2 = 0;
    Eigen::MatrixXd V_box(1,3);
    ifstream input;
    input.open(fileCgaL);
    while (true) {
        float x1,y1,z1;
        input>>x1>>y1>>z1;
        if( input.eof() ) break;
        push_backMat(V_box, Vector3d(x1,y1,z1), cont2);
        cont2++;
    }

    input.close();
    return V_box;
}
vector <vector<float>> Process::readEigensLap(){
    string typeFileMatrix;
    switch(typeAdjMatrix){
        case 'L': typeFileMatrix = "_eigenLap.cgal";
            break;
        case 'C': typeFileMatrix = "_eigen.cgal";
            break;
        case 'N': typeFileMatrix = "_eigenCorrDist.cgal";//"_eigenCorr.cgal";
            break;
        case 'A': typeFileMatrix = "_eigenWeight.cgal";
            break;
        case 'W': typeFileMatrix = "_eigenWeight_W.cgal";
            break;
        default :
             cout << "Invalid type of matrix read eigen" << endl;
    }
    cout << typeFileMatrix << endl;
    cout << "changeNumberVert " << changeNumberVert << endl;
    std::string dirNewDesc;
    if(changeNumberVert == true){
        dirNewDesc = PATH_MESHES_400;
    }
    else {
        if(typeAdjMatrix == 'N'){
            cout << "entro" <<endl;
            dirNewDesc = "/mnt/data/Code/KeyPoints/project/SampleMeshes/db/newDB_Skel/T";
        }
        else{
            dirNewDesc = PATH_MESHES;
        }
    }
    cout <<typeAdjMatrix<<endl;
    cout << dirNewDesc << endl;
    vector <vector<float>> descData;
    for (int i = 0; i < numberShapes; i++){
        string name = "T" + to_string(i);
        std::string dirDesc = dirNewDesc + string(to_string(i)) + string("/") + name + typeFileMatrix;
//        if(changeNumberVert == true){
//            dirDesc = PATH_MESHES_400 + string(to_string(i)) + string("/") + name + typeFileMatrix;
//        }
//        else{
//            dirDesc = PATH_MESHES + string(to_string(i)) + string("/") + name + typeFileMatrix;
//        }
      //  cout << dirDesc << endl;
      ///////////////////////////////////  const char * fileDirDesc = dirDesc.c_str();
        string fileDirDesc = dirDesc.c_str();
        vector<float> desc;
//        cout << desc.size() << endl;
        {
            std::ifstream in(fileDirDesc);
            while (true) {
                float value;
                in >> value;
                if( in.eof() ) break;
                desc.push_back(value);
            }
        }
        descData.push_back(desc);
//        cout << descData.size() << endl;
      //  cout << descData[i].size() << endl;
    }
    return descData;
}
Eigen::MatrixXd Process::readVertForGraph(){
    string fileVertices;
    string  fileEdges;
    if(typeSkeleton == true){
        fileVertices = "vertex.cgal";
        fileEdges = "skel2.cgal";
    }
    else {
        fileVertices = "vertex-sm.cgal";
        fileEdges = "skel-sm.cgal";
    }
    Eigen::MatrixXd vertices = readEdges(fileVertices);
    Eigen::MatrixXd edgesConect = readEdges(fileEdges);
    vector<int> edgesItem(edgesConect.rows(),-1);
    for (int i = 0; i < edgesConect.rows() ; i ++){
        for (int j = 0; j < vertices.rows() ; j ++){
            if (edgesConect.row(i) == vertices.row(j)){
                edgesItem[i]= j;
                break;
            }
        }
    }
    Eigen::MatrixXd graph = getLapGraph(vertices,edgesItem);
    return graph;
}
vector<float> Process::distSurfVert(Eigen::MatrixXd fileVertices, string filevertCorrSurf){
    cout << currDiag <<endl;
    Eigen::MatrixXd vertices = fileVertices;
    Eigen::MatrixXd vertCorrSurf = readEdges(filevertCorrSurf);//"correspondance-poly.cgal");
    vector<float> numVertSurf(vertices.rows(),0);
    int n=0;
    float tmpDist;
    for (int i = 0; i < vertices.rows(); i ++){
        int count = 0;
        float MinDist = std::numeric_limits<float>::infinity();
        float MaxDist = 0.0;
        for (int j = n; j < vertCorrSurf.rows() ; j= j+2){
          //  cout << vertices.row(i) << vertCorrSurf.row(j)<<endl;
            if (vertices.row(i) == vertCorrSurf.row(j)){
                count+=1;
                float newDist = sqrt((vertices(i,0) - vertCorrSurf(j+1,0))*(vertices(i,0) - vertCorrSurf(j+1,0))+
                (vertices(i,1) - vertCorrSurf(j+1,1))*(vertices(i,1) - vertCorrSurf(j+1,1))+
                (vertices(i,2) - vertCorrSurf(j+1,2))*(vertices(i,2) - vertCorrSurf(j+1,2)));
//                cout<< newDist<<endl;
//                cout << vertices(i,0) << " "<< vertices(i,1) << " "<<vertices(i,2) <<endl;
//                cout << vertCorrSurf(j+1,0) << " "<< vertCorrSurf(j+1,1) << " "<<vertCorrSurf(j+1,2) <<endl;
                if(MinDist > newDist){
                    MinDist = newDist;
                }
                if(MaxDist < newDist){
                    MaxDist = newDist;
                }
                tmpDist = MinDist/MaxDist;
                //cout << "entro"<<endl;
            }
            else{ n += (count *2);
               //  cout << "min: "<<MinDist<<endl;
               //  cout << "max: "<<MaxDist<<endl;
                float dist =  MinDist/MaxDist;
                if(count == 0 || isinf(dist) || isnan(dist)){
                    numVertSurf[i]= tmpDist; break;
                }
                else{
                 numVertSurf[i]= dist; break;}
            }

        }
        numVertSurf[i]= tmpDist;
    }
//    vector<float>::iterator it = numVertSurf.begin();
//    int i = 1;
//    for (; it != numVertSurf.end(); it++ ){
//       cout<< i <<":" <<*it<< endl;
//       i++;
//    }
    return numVertSurf;
}
vector<int> Process::countSurfVert(Eigen::MatrixXd fileVertices, string filevertCorrSurf){
    cout << currDiag <<endl;
    Eigen::MatrixXd vertices = fileVertices;
    Eigen::MatrixXd vertCorrSurf = readEdges(filevertCorrSurf);//"correspondance-poly.cgal");
    vector<int> numVertSurf(vertices.rows(),0);
    int n=0;
    for (int i = 0; i < vertices.rows(); i ++){
        int count = 0;
        for (int j = n; j < vertCorrSurf.rows() ; j= j+2){
          //  cout << vertices.row(i) << vertCorrSurf.row(j)<<endl;
            if (vertices.row(i) == vertCorrSurf.row(j)){
                count+=1;
                numVertSurf[i]= count;
               // cout << "entro"<<endl;
            }
            else{ n += (count *2); break;}
        }
    }
//    vector<int>::iterator it = numVertSurf.begin();
//    int i = 1;
//    for (; it != numVertSurf.end(); it++ ){
//       cout<< i <<":" <<*it<< endl;
//       i++;
//    }
    return numVertSurf;
}

vector<float> Process::getActLapEigen(){
    Eigen::MatrixXd graph = readVertForGraph();
    EigenSolver<MatrixXd> es(graph);
    vector<float> E;
    for(int i = 0; i < graph.rows(); ++i){
        E.push_back(es.eigenvalues().col(0)[i].real());
    }
    sort(E.begin(),E.end());
    /////////////////////const char *nameEigen = getEigenFileName(currentShape);
    string nameEigen = getEigenFileName(currentShape);
    std::ofstream outputLapEigen(nameEigen);
    for (int i =0; i< E.size(); i++){
        outputLapEigen <<E[i]<< "\n";
    }
    outputLapEigen.close();
    return E;
}
map<int,int> Process::getJunction(Eigen::MatrixXd graph){
    map <int, int> junct;
    for (int i = 0; i < graph.rows(); i++){
        int n = 0;
        for(int j = 0; j < graph.rows(); j++){
            if(graph(i,j) == 1 )
            {
                n= n +1;
            }
        }
        if(junct.count(n) > 0){
            junct[n] = junct[n] +1;
        }
        else{
            junct[n] = 1;
        }
    }
    return junct;
}
map<int,int> Process::getCurJunct(){
    Eigen::MatrixXd graph = readVertForGraph();
    map<int, int> descJunct = getJunction(graph);
    string chosenFileJunct = getResultFileName(currentShape) + "_junct.cgal";
    ///////////////const char *nameJunct = chosenFileJunct.c_str();
    string nameJunct = chosenFileJunct.c_str();
    std::ofstream outputJunct(nameJunct);
    map<int,int>::iterator it = descJunct.begin();
    for (; it != descJunct.end(); it++ ){
       outputJunct << it->first <<" " << it -> second<< endl;
    }
    outputJunct.close();
    return descJunct;
}
vector <map<int,int>> Process::readJunct(){
    vector <map<int,int>>  descData;
    for (int i = 0; i < numberShapes; i++){
        string name = "T"+string(to_string(i));
        string dirDesc = PATH_MESHES + string(to_string(i))+string("/")+name +"_junct.cgal";
        ///////////////////const char * fileDirDesc = dirDesc.c_str();
        string fileDirDesc = dirDesc.c_str();
        map<int,int> desc;
        {
            std::ifstream in(fileDirDesc);
            while (true) {
                int key, value;
                in >> key >>value;
                if( in.eof() ) break;
                desc[key] = value;
            }
        }
        descData.push_back(desc);
    }
    return descData;
}
vector<map<int,int>> Process::compareJunct(map<int,int> actJunct, vector<map<int,int>> dataJunct){
    vector<map<int,int>> dif;
    for(int i = 0; i < dataJunct.size(); i++){
        map<int,int>::iterator itData;
        map<int,int>::iterator itAct = actJunct.begin();
        map<int,int>comp;
        for (; itAct != actJunct.end(); itAct++ ){
            itData = dataJunct[i].find(itAct->first);
            if(itData != dataJunct[i].end() && itAct->first != 2 )
                comp[itData->first] = abs(itData->second - itAct->second);
            else if(itData == dataJunct[i].end() && itAct->first != 2 ){
                comp[itAct->first] = itAct->second;
            }
        }
        for(itData = dataJunct[i].begin(); itData != dataJunct[i].end(); itData++){
            itAct = actJunct.find(itData->first);
            if(itAct == actJunct.end()){
                    comp[itData->first] = itData->second;
            }
        }
        dif.push_back(comp);
    }
    return dif;
}
vector<pair<float, string>> Process::sumJunct(vector<map<int,int>> dataJunct){
    map<string,float> result;
    vector<pair<float, string>> mapVector;
    for (int i = 0; i < dataJunct.size(); i ++){
        map<int,int>::iterator itData = dataJunct[i].begin();
        int sumData=0;
        for (; itData != dataJunct[i].end(); itData++){
                sumData += itData->second;
        }
        pair<float, string> p1;
        p1.first = sumData;
        p1.second = "T"+to_string(i);//shapes[i];
        mapVector.push_back(p1);
        result["T"+to_string(i)] = sumData;//shapes[i]
    }
    sort(mapVector.begin(), mapVector.end());
    return mapVector;

}
map<string, float> Process::compareLargerDesc(vector<vector<float>> descData, vector<float> descAct){
    int actSize = descAct.size();
    map<string,float> result;
    for (int i = 0; i < descData.size(); i ++){
        int actDataSize = descData[i].size();
        int descSize;
        if (actSize < numberEigen || actDataSize < numberEigen){
            if(actSize <= actDataSize){
                descSize = actSize;
            }
            else{
                descSize = actDataSize;
            }
        }
        else{
            descSize = numberEigen;

        }
       // cout << "numberEigen: "<< descSize<<endl;
        vector<float> v1(descData[i].end()- descSize, descData[i].end() -1);
        //for (int i = 0; i < v1.size(); i++){
          //cout << "v1: " << v1[i]<<endl;
        //}
        //cout << "--------------------------------------" << endl;
        vector<float> v2(descAct.end()- descSize, descAct.end() -1);
        //result["T" + to_string(i)] = vectors_distance(v2,v1);
        result[to_string(i)] = vectors_distance(v2,v1);
    }

    return result;
}
map<string, float> Process::compareSmallerDesc(vector<vector<float>> descData, vector<float> descAct){
    //cout << "entro smaller descriptor"<<endl;
    int actSize = descAct.size();
    map<string,float> result;
    for (int i = 0; i < descData.size(); i ++){
        int actDataSize = descData[i].size();
        int descSize;
        if (actSize < numberEigen || actDataSize < numberEigen){
            if(actSize <= actDataSize){
                descSize = actSize;
            }
            else{
                descSize = actDataSize;
            }
        }
        else{
            descSize = numberEigen;

        }
       // cout << "numberEigen: "<< descSize<<endl;
        vector<float> v1(descData[i].begin(), descData[i].begin() + descSize);
        // for (int i = 0; i < v1.size(); i++){
        //   cout << "v1: " << v1[i]<<endl;
        // }
        // cout << "--------------------------------------" << endl;
        vector<float> v2(descAct.begin(), descAct.begin() + descSize);
        //result["T" + to_string(i)] = vectors_distance(v2,v1);
        result[to_string(i)] = vectors_distance(v2,v1);
    }
    return result;
}
int Process::startDescript(vector<float> descAct){
  int startDescAct = 0;
  while(descAct[startDescAct] < 0){
    startDescAct= startDescAct +1;
  }
  return startDescAct;
}
map<string, float> Process::compareSmallerPosDesc(vector<vector<float>> descData, vector<float> descAct){
    //cout << "entro smaller descriptor"<<endl;
    int actSize = descAct.size();
    map<string,float> result;

    int startDescAct = startDescript(descAct);


    for (int i = 0; i < descData.size(); i ++){
        int actDataSize = descData[i].size();
        int descSize;
        if (actSize < numberEigen || actDataSize < numberEigen){
            if(actSize <= actDataSize){
                descSize = actSize;
            }
            else{
                descSize = actDataSize;
            }
        }
        else{
            descSize = numberEigen;

        }
       // cout << "numberEigen: "<< descSize<<endl;
        int startDescData = startDescript(descData[i]);
        vector<float> v1(descData[i].begin() + startDescData, descData[i].begin() + startDescData + descSize);
        // for (int i = 0; i < v1.size(); i++){
        //   cout << "v1: " << v1[i]<<endl;
        // }
        // cout << "--------------------------------------" << endl;
        vector<float> v2(descAct.begin() + startDescAct, descAct.begin() + startDescAct+ descSize);
        //result["T" + to_string(i)] = vectors_distance(v2,v1);
        result[to_string(i)] = vectors_distance(v2,v1);
    }
    return result;
}
void Process::compEucDist(){
    vector <vector<float>> matDesc = readEigensLap();
    vector<float> actDes = getActLapEigen();
    map<string, float> eucDist = compareLargerDesc(matDesc, actDes);
    map<string,float>::iterator it = eucDist.begin();
    vector<pair<float, string>> mapVector;
    for (; it != eucDist.end(); it++ ){
        pair<float, string> p1;
        p1.first=it -> second;
        p1.second=it -> first;
        mapVector.push_back(p1);
    }
    sort(mapVector.begin(), mapVector.end());
    vector<pair<float,string>>::iterator itS = mapVector.begin();
    string result = getResultFileName(currentShape)+".res";
    const char *nameResult = result.c_str();
    std::ofstream outputLapResult(nameResult);
    for (; itS != mapVector.end(); itS++ ){
          cout << itS->first <<": " << itS -> second<< endl;
          outputLapResult <<itS->first <<": " << itS -> second<< "\n";
    }
    outputLapResult.close();
}

void Process::compMatrixDissim(){
    vector <vector<float>> matDesc = readEigensLap();
    cout << "ok"<<endl;
    string result = getResultFileName(currentShape)+"_"+typeAdjMatrix+"_MatrixDist.res";
    cout << result<< endl;
    ///////const char *matrixResult = result.c_str();
    string matrixResult = result.c_str();
    std::ofstream outputMatResult(matrixResult);
    for(int i = 0; i < matDesc.size(); i++){
        vector<float> actDes = matDesc[i];
        //cout << actDes.size() << endl;
        map<string, float> eucDist = compareSmallerPosDesc(matDesc, actDes);
      //  map<string, float> eucDist = compareLargerDesc(matDesc, actDes);
        //map<string, float> eucDist = compareSmallerDesc(matDesc, actDes);

        vector<pair<int, float>> mapVector = sortMap(eucDist);

        vector<pair<int,float>>::iterator it = mapVector.begin();
//        for (; it != eucDist.end(); it++ ){
//            outputMatResult << it -> second << " ";
//        }
        for (; it != mapVector.end(); it++ ){
            outputMatResult << it -> second << " ";
        }
        outputMatResult << "\n";
    }
    outputMatResult.close();
}

vector<pair<int, float>> Process::sortMap(map<string, float> eucDist){
    map<string,float>::iterator it = eucDist.begin();
    vector<pair<int, float>> mapVector;
    for (; it != eucDist.end(); it++ ){
        pair<int, float> p1;
        p1.first= std::stoi (it -> first,nullptr);
        p1.second=it -> second;
        mapVector.push_back(p1);
    }
    sort(mapVector.begin(), mapVector.end());
    return mapVector;
}
void Process::compJunct(){
    map<int,int> descJunct = getCurJunct();
    vector<map<int,int>> dataJunct= readJunct();
    vector<map<int,int>> compare = compareJunct(descJunct, dataJunct);
//    cout << "print compare: ----------------------------" << endl;
//    for(int ii = 0; ii < compare.size(); ii++){
//        cout<<"modelos "<< to_string(ii) <<endl;
//        map<int,int>::iterator itCom = compare[ii].begin();
//        for (; itCom != compare[ii].end(); itCom++ ){
//            cout << itCom->first <<": " << itCom -> second<< endl;
//        }
//     }
    cout << "Result: ----------------------------" << endl;
    vector<pair<float, string>> compareDist = sumJunct(compare);
    vector<pair<float, string>>::iterator itDist = compareDist.begin();
    for (; itDist != compareDist.end(); itDist++ ){
        cout << itDist->first <<": " << itDist -> second<< endl;
    }
}

void Process::saveGraph(Eigen::MatrixXd graph){
    string nameGraph = getShapeName(currentShape);
    std::ofstream outputLapEigen(nameGraph);
}
void Process::printActualEigenLaplace(vector<float> actDes){
    cout << "--------------Actual Eigen Laplace--------------"<<endl;
    for (int i =0; i< actDes.size(); i++){
        cout<<actDes[i]<< endl;
    }
}
void Process::printDataBaseEigenLaplace(){
    vector <vector<float>> matDesc = readEigensLap();
    cout << "--------------Data base Eigen Laplace--------------"<<endl;
    for (int i =0; i< matDesc.size(); i++){
        cout << "--------Model: " << i << "--------" <<endl;
        for(int j= 0; j <matDesc[i].size();j++){
            cout<<matDesc[i][j]<< endl;
        }
    }
}
void Process::printCurrentJunct(){
    map<int,int> descJunct = getCurJunct();
    map<int,int>::iterator it = descJunct.begin();
    for (; it != descJunct.end(); it++ ){
        cout << it->first <<": " << it -> second<< endl;
    }
}
void Process::printDataBaseJunct(){
    vector<map<int,int>> dataJunct= readJunct();
    cout << "data: ----------------------------" << endl;
    for(int ii = 0; ii < dataJunct.size(); ii++){
        cout<<"modelos T"<< to_string(ii)  <<endl; //shapes[ii]
        map<int,int>::iterator itCom = dataJunct[ii].begin();
        for (; itCom != dataJunct[ii].end(); itCom++ ){
            cout << itCom->first <<": " << itCom -> second<< endl;
        }
     }
}


void Process::saveEigen (string chosenFile, char adjMatrix){
    typeAdjMatrix = adjMatrix;
    string chosenFileV = chosenFile + "_vertex.cgal";
    const char * fileVertices = chosenFileV.c_str();
    Eigen::MatrixXd vertices = readEdges(fileVertices);
    string chosenFileE = chosenFile + "_skel.cgal";
    const char  *fileEdges = chosenFileE.c_str();
    Eigen::MatrixXd edgesConect = readEdges(fileEdges);
    cout << "edgesConect.rows(): "<<edgesConect.rows()<< endl;
    vector<int> edgesItem(edgesConect.rows(),-1);
    for (int i = 0; i < edgesConect.rows() ; i ++){
        for (int j = 0; j < vertices.rows() ; j ++){
            if (edgesConect.row(i) == vertices.row(j)){
                edgesItem[i] = j;
                //cout << "edgesItem: "<<edgesItem[i]<< endl;
                break;
            }
        }
    }
    Eigen::MatrixXd graph(vertices.rows(), vertices.rows());
    currDiag = chosenFile + "_correspondance-poly.cgal";
    cout << currDiag<< endl;
    graph = getLapGraph(vertices,edgesItem);

    //cout << graph << endl;
    EigenSolver<MatrixXd> es(graph);
    vector<float> E;
    vector<float> V;
    float sum = 0.0;
    for(int i = 0; i < graph.rows(); ++i){
        E.push_back(es.eigenvalues().col(0)[i].real());
        for (int j = 0; j < es.eigenvectors().col(i).rows() ; j ++){
            sum = sum + es.eigenvectors().col(i)[j].real();
        }
        V.push_back(sum);
        sum=0.0;
    }
    map<float,float> featureVector;
    float eigenValue;
    for (int i = 0; i < E.size(); i++){
        //eigenValue = abs(E[i]);
        eigenValue = E[i];
        if (featureVector.find(eigenValue) != featureVector.end()){
            featureVector[eigenValue] += V[i];
        }
        else{
            featureVector[eigenValue] = V[i];
       }
    }
    sort(E.begin(),E.end());
    string typeFileSaveMatrix;
    switch(typeAdjMatrix){
        case 'L': typeFileSaveMatrix = "_eigenLap.cgal";
            break;
        case 'C': typeFileSaveMatrix = "_eigen.cgal";
            break;
        case 'N': //typeFileSaveMatrix = "_eigenCorr.cgal";
                typeFileSaveMatrix = "_eigenCorrDist.cgal";
            break;
        case 'A': //typeFileSaveMatrix = "_eigenCorr.cgal";
                typeFileSaveMatrix = "_eigenWeight.cgal";
            break;
        case 'W': //typeFileSaveMatrix = "_eigenCorr.cgal";
                typeFileSaveMatrix = "_eigenWeight_W.cgal";
            break;
        default :
             cout << "Invalid type of matrix" << endl;
    }
    string chosenFileEig = chosenFile + typeFileSaveMatrix;
    const char *nameEigen = chosenFileEig.c_str();
    std::ofstream output2(nameEigen);
    for (int i = 0; i <E.size(); i ++){
          output2 << E[i] << "\n";

          //cout << "E: " << eigenValue << endl;
    }
    output2.close();
    //saveGraph
    string chosenFileG = chosenFile + "_graph_"+ typeFileSaveMatrix;
    const char *nameG = chosenFileG.c_str();
    std::ofstream outputG(nameG);
    outputG << graph << "\n";

    //Save descriptor as GISIFs
//    string chosenFileDesc = chosenFile + typeFileSaveMatrix +"_desc.cgal";
//    const char *nameDesc = chosenFileDesc.c_str();
//    std::ofstream outputDesc(nameDesc);
//    for(map<float,float>::iterator it = featureVector.begin(); it != featureVector.end(); ++it){
//          outputDesc << it->second << "\n";
//    }
//  //  outputDesc << std::fixed << std::setprecision(8) << value;
//    outputDesc.close();
}

//void saveJunct(string chosenFile){
//    string chosenFileV = chosenFile + "_vertex.cgal";
//    const char * fileVertices = chosenFileV.c_str();
//    Eigen::MatrixXd vertices = readEdges(fileVertices);
//    string chosenFileE = chosenFile + "_skel.cgal";
//    const char  *fileEdges = chosenFileE.c_str();
//    Eigen::MatrixXd edgesConect = readEdges(fileEdges);
//    vector<int> edgesItem(edgesConect.rows(),-1);
//    for (int i = 0; i < edgesConect.rows() ; i ++){
//        for (int j = 0; j < vertices.rows() ; j ++){
//            if (edgesConect.row(i) == vertices.row(j)){
//                edgesItem[i]= j;
//                break;
//            }
//        }
//    }
//    Eigen::MatrixXd graph = getGraph(vertices,edgesItem);
//    map<int, int> descJunct = getJunction(graph);

//    string chosenFileJunct = chosenFile + "_junct.cgal";
//    const char *nameJunct = chosenFileJunct.c_str();
//    std::ofstream output(nameJunct);
//    map<int,int>::iterator it = descJunct.begin();
//    for (; it != descJunct.end(); it++ ){
//       output << it->first <<" " << it -> second<< endl;
//    }
//    output.close();
//}
