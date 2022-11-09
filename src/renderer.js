const { ipcRenderer } = require("electron");
const { sqrt } = require("mathjs");
const math = require("mathjs");

function run() {
let n = document.getElementById('n').value;
let sum_Xi = document.getElementById('sum_Xi').value;
let sum_Zi = document.getElementById('sum_Zi').value;
let sum_Yi = document.getElementById('sum_Yi').value;
let sum_Xi_sqr = document.getElementById('sum_Xi_sqr').value;
let sum_Zi_sqr = document.getElementById('sum_Zi_sqr').value;
let sum_Yi_sqr = document.getElementById('sum_Yi_sqr').value;
let sum_XiZi = document.getElementById('sum_XiZi').value;
let sum_YiZi = document.getElementById('sum_YiZi').value;
let sum_XiYi = document.getElementById('sum_XiYi').value;
let k = document.getElementById('k').value;
let Xo = document.getElementById('Xo').value;
let Zo = document.getElementById('Zo').value;
let hằng_số_t_alpha = document.getElementById('hằng_số_t_alpha').value;

const matrix1 = [[n, sum_Xi ,sum_Zi ], [sum_Xi ,sum_Xi_sqr ,sum_XiZi ], [sum_Zi ,sum_XiZi , sum_Zi_sqr ]]
const matrix_XtX = math.matrix(matrix1)    

const matrix2 = [sum_Yi , sum_XiYi , sum_YiZi ]
const matrix_XtY = math.matrix(matrix2)

const det = math.det(matrix_XtX)

const inverse = (m) => {
    return math.inv(m);
  };
  const matrix_XtX_inverse = inverse(matrix_XtX);

  const betaMu = math.multiply(matrix_XtX_inverse, matrix_XtY)

  document.getElementById('result1').innerHTML = betaMu

let sum_resid_sqr = sum_Yi_sqr - (betaMu.subset(math.index(0))*sum_Yi + betaMu.subset(math.index(1))*sum_XiYi + betaMu.subset(math.index(2))*sum_YiZi);
let sigmaMu_sqr = sum_resid_sqr/(n-k);
const khoảng_tin_cậy_betaMu = []
const standardError_array = []

for (j = 1; j < betaMu.length; j++){
    let variance_betaMu = sigmaMu_sqr*(matrix_XtX_inverse[j][j])
    let standardError = Math.sqrt(variance_betaMu)
    standardError_array.push(standardError)
    let KTC_left = betaMu - hằng_số_t_alpha*standardError
    let KTC_right = betaMu + hằng_số_t_alpha*standardError
    khoảng_tin_cậy_betaMu.push(KTC_left + '-' + KTC_right)
}
document.getElementById('result2').innerHTML = sigmaMu_sqr

const matrix_Xto = math.matrix([1, Xo, Zo])
const matrix_Xo = math.matrix([[1], [Xo], [Zo]])
let Yo = betaMu.subset(math.index(0)) + betaMu.subset(math.index(1))*Xo + betaMu.subset(math.index(2))*Zo
let variance_YoMu = sigmaMu_sqr*math.multiply(math.multiply(matrix_Xto, matrix_XtX_inverse), matrix_Xo)
let standardError_YoMu = Math.sqrt(variance_YoMu)
let KTC_left_Yo = Yo - hằng_số_t_alpha*standardError_YoMu;
let KTC_right_Yo = Yo + hằng_số_t_alpha*standardError_YoMu;
const KTC_Y_Xo = [KTC_left_Yo, KTC_right_Yo]

document.getElementById('result3').innerHTML = Yo 
}
