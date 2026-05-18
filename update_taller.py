import re

file_path = "/home/juan-calixto/Documents/programming/personal/universidad/gerencia_y_gestion/informe/taller1.tex"

with open(file_path, "r", encoding="utf-8") as f:
    lines = f.readlines()

out = []
found = False
for line in lines:
    if line.startswith(r"\doctype{Documentation Template}"):
        found = True
        out.append(r"""\doctype{Informe de Gestión}
\title{Análisis de Buenas y Malas Prácticas en Proyectos de Ingeniería: Estudio de Casos SECOP}

%----------------------------------------------------------
% Authors, Affiliations and dates
%----------------------------------------------------------

\author[a,1]{Juan Felipe Calixto Londoño}


%----------------------------------------------------------

\affil[a]{Universidad Nacional de Colombia, Facultad de Ingeniería}

%----------------------------------------------------------

\dates{Semestre 2026-1}

%----------------------------------------------------------
% Corresponding author-, Document- information
%----------------------------------------------------------

\corres{\textsuperscript{\textasteriskcentered}Autores correspondientes: 
\href{mailto:jcalixtol@unal.edu.co}{jcalixtol@unal.edu.co}


\journalname{LATEX}
\journal{Rho-class (v.3) 2026}
\theday{\today}
\vol{X}
\no{X}

%----------------------------------------------------------
% Abstract and Keywords
%----------------------------------------------------------

\begin{abstract}
    Este documento presenta un análisis crítico sobre las buenas y malas prácticas en la formulación, ejecución y seguimiento de proyectos de ingeniería en el marco de la contratación pública en Colombia. Se utilizan como base de estudio de caso dos contratos extraídos del Sistema Electrónico de Contratación Pública (SECOP). Para ilustrar las buenas prácticas, se examina el contrato SA-MEN-05-2022, destacando su transparencia, planificación y cumplimiento. Por otro lado, se analiza el contrato IDU-1635-2019 como ejemplo de malas prácticas, evidenciando problemas de sobrecostos, retrasos y deficiencias en la planeación. El propósito de este informe es establecer lecciones aprendidas aplicables a la gerencia de proyectos de ingeniería.
\end{abstract}

%----------------------------------------------------------

\keywords{SECOP, Gerencia de Proyectos, Contratación Pública, Buenas Prácticas, Ingeniería}

%----------------------------------------------------------

\begin{document}

    \maketitle
    \thispagestyle{firststyle}

%----------------------------------------------------------

\section{Introducción}

    La gerencia de proyectos de ingeniería en el sector público requiere un control estricto, transparencia y una planeación rigurosa para garantizar el uso adecuado de los recursos del Estado. En Colombia, el Sistema Electrónico de Contratación Pública (SECOP) es la plataforma principal que registra y permite la trazabilidad de estos contratos, convirtiéndose en una fuente invaluable de información para analizar el desempeño de la gestión pública.
    
    En este documento se contrastan dos proyectos de infraestructura y servicios de ingeniería. El primero, el contrato \textbf{SA-MEN-05-2022}, representa un estándar de buenas prácticas, mientras que el contrato \textbf{IDU-1635-2019} sirve como caso de estudio de las deficiencias, retrasos y sobrecostos comunes frente a la falta de rigor técnico y gerencial abordado durante las etapas tempranas del ciclo de vida del proyecto.

\section{Marco de Referencia: Contratación Estatal y SECOP}

    El SECOP incentiva la transparencia administrativa. Para evaluar el éxito o fracaso de los contratos, se utilizan criterios propios de la gestión de proyectos, tales como:
    \begin{itemize}
        \item \textbf{Alcance:} Definición clara de los requerimientos desde los pliegos de condiciones.
        \item \textbf{Cronograma y Tiempo:} Cumplimiento estricto de los hitos y del plazo de ejecución final estipulado.
        \item \textbf{Costo:} Adherencia al presupuesto oficial, reduciendo el número de adiciones financieras ("otrosí").
        \item \textbf{Gestión de Riesgos:} Identificación temprana de amenazas técnicas, ambientales o sociales.
    \end{itemize}

\section{Análisis de Buenas Prácticas: Contrato SA-MEN-05-2022}

    El contrato \textbf{SA-MEN-05-2022} (Ministerio de Educación Nacional) se destaca como un referente positivo de gestión. Del análisis de los documentos en SECOP, se extraen las siguientes buenas prácticas:

    \subsection{Planificación y Estudios Previos Rigurosos}
    La documentación resalta una solidez notable en los estudios previos. El análisis exhaustivo de mercado permitió estructurar un presupuesto realista, definiendo requerimientos técnicos alcanzables y disminuyendo la incertidumbre antes de abrir la licitación, asegurando proyectos de infraestructura educativa viables.

    \subsection{Gestión de Riesgos y Pólizas}
    La matriz de riesgos del proyecto estuvo bien estructurada, con una distribución equitativa de las cargas entre la entidad contratante y el contratista. Se dispusieron pólizas adecuadas para cubrir aspectos esenciales de cumplimiento, pago de salarios y estabilidad de la obra civil.

    \subsection{Transparencia y Criterios Claros de Adjudicación}
    Los pliegos de condiciones se diseñaron garantizando la pluralidad de oferentes mediante puntajes sensatos, evitando "pliegos sastre" o direccionamientos contractuales indebidos.

\section{Análisis de Malas Prácticas: Contrato IDU-1635-2019}

    En contraste, el contrato \textbf{IDU-1635-2019}, centrado en el desarrollo de obra pública en Bogotá (Instituto de Desarrollo Urbano), ejemplifica una cadena de anomalías y malas prácticas de contratación:

    \subsection{Deficiencias en los Diseños y Maduración del Proyecto}
    El principal detonante de los fracasos de este contrato radicó en iniciar la fase de ejecución sin diseños estructurados ni predios liberados y saneados. Comenzar la etapa constructiva con incertidumbres masivas conllevó obligatoriamente a constantes interrupciones, paralizaciones y suspensiones temporales en las diferentes zonas de obra.

    \subsection{Sobrecostos y Modificaciones Constantes (Otrosí)}
    Debido a una planeación deficiente, como el descubrimiento tardío de redes de servicios públicos no inventariadas, el proyecto sufrió múltiples rediseños de ingeniería durante la ejecución. Como resultado, se aprobaron prórrogas desmedidas y varias adiciones presupuestales de gran tamaño que incrementaron notoriamente el pago a favor del contratista.

    \subsection{Debilidades en la Interventoría y Seguimiento}
    En auditorías se ha demostrado la falta de articulación entre los consorcios de obra civil y la interventoría para ejecutar planes de contingencia oportunos y recuperar los retrasos en la ruta crítica del cronograma planificado.

\section{Conclusiones}

    El análisis comparativo destaca que la causa fundamental del éxito en el contrato \textbf{SA-MEN-05-2022} descansó en la madurez y exhaustividad de su fase de iniciación/planeación. Al contar con un alcance sólidamente delimitado y riesgos técnicos correctamente ponderados, la ejecución transcurrió dentro de los límites del control gerencial.
    
    Por el contrario, el contrato \textbf{IDU-1635-2019} evidencia las repercusiones financieras y sociales de acelerar fases de preinversión licitando obras con diseños esquemáticos incompletos y con externalidades no dirimidas por el estado. Estas fallas comunes no solo castigan al erario a través de billonarios sobrecostos, sino que terminan afectando drásticamente el bienestar de la ciudadanía.

\end{document}
""")
        break
    else:
        out.append(line)

with open(file_path, "w", encoding="utf-8") as f:
    f.writelines(out)

